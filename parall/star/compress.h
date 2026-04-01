#ifndef COMPRESS_H
#define COMPRESS_H

#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
extern int errno;
#include <list>
#include <string>


class ckpt_compressor {
private:
  enum { MAX_FILE_PATH = 1024 };
  pid_t cpid;
  std::string pattern;
  int freq;
  int rfd, wfd;

  int add_to_tar(char *tarname, char *ckptname){
    char command[2*MAX_FILE_PATH+10]; // approximately MAX size of the command
    sprintf(command,"tar -rvf %s %s 2>/dev/null 1>/dev/null",tarname, ckptname);
    int status = system(command);
    if(  WEXITSTATUS(status) != 0 ){
      printf("ckpt_compressor: Error adding ckpt %s into the tar %s. Keep ckpt file.", ckptname, tarname);
      return -1;
    }

    if( unlink(ckptname) ){
      char buf[MAX_FILE_PATH + 100];
      sprintf(buf,"ckpt_compressor: unlink(%s)",ckptname);
      perror(buf);
      return -1;
    }

    return 0;
  }

  int gzip_tar(char *tarname){
    char command[2*MAX_FILE_PATH+10]; // approximately MAX size of the command
    sprintf(command,"gzip -fq %s 1>/dev/null 2>/dev/null",tarname);
    int status = system(command);
    if(  WEXITSTATUS(status) != 0 ){
      printf("ckpt_compressor: Error compressing tar %s.", tarname);
      return -1;
    }

    return 0;
  }


  void compressor_thread(){
    std::list<std::string> l;
    char buf[MAX_FILE_PATH];
    int ckpt_cnt = 0, tar_cnt = 0;
    int cnt;
    char tar_name[MAX_FILE_PATH];
    sprintf(tar_name, "%s_%d.tar",pattern.c_str(), tar_cnt);
    while( (cnt = read(rfd, buf, MAX_FILE_PATH) ) > 0 ){
      char *ptr = buf;
      for(int i = 0; i < cnt; i++){
        if(buf[i] == '\0' ){
          if( add_to_tar(tar_name, ptr) ){
            printf("ckpt_compressor: Cannot add ckpt %s to the tar %s. Skip.", ptr, tar_name);
          }
          //printf("Read ckpt %s, count = %d, freq = %d\n", ptr, ckpt_cnt, freq);
          ptr = buf + i + 1;
          i++;
          ckpt_cnt++;
          if( ckpt_cnt >= freq ){
            printf("Compress tar %s\n", tar_name);
            gzip_tar(tar_name);
            // Switch to the next tar
            tar_cnt++;
            ckpt_cnt = 0;
            sprintf(tar_name, "%s_%d.tar",pattern.c_str(), tar_cnt);
            // Cleanup stalled files if any
            unlink(tar_name);
            char tmp_buf[MAX_FILE_PATH];
            sprintf(tmp_buf,"%s.gz",tar_name);
            if( unlink(tar_name) ){
              if( errno != ENOENT ){
                char buf[MAX_FILE_PATH + 100];
                sprintf(buf,"ckpt_compressor: unlink(%s)",tar_name);
                perror(buf);
              }
            }
            if( unlink(tmp_buf) ){
              if( errno != ENOENT ){
                char buf[MAX_FILE_PATH + 100];
                sprintf(buf,"ckpt_compressor: unlink(%s)",tmp_buf);
                perror(buf);
              }
            }

          }
        }
      }
    }

    if( ckpt_cnt ){
      // Compress remaining ckpts
      gzip_tar(tar_name);
    }
  }

  bool state_ok;
public:
  ckpt_compressor(){
    cpid = -1;
    pattern = "";
    state_ok = false;
  }

  int init(const char *_pattern, int _freq)
  {
    pattern = _pattern;
    freq = _freq;
    int chan[2];

    if( pipe(chan) ){
      perror("Compressor pipe");
      state_ok = false;
      return -1;
    }
    rfd = chan[0];
    wfd = chan[1];

    cpid = fork();
    if( cpid < 0 ){
      perror("Compressor fork");
      state_ok = false;
      return -1;
    }
    if( cpid == 0 ){
      close(wfd);
      compressor_thread();
      // Exit after we done
      exit(0);
    }
    close(rfd);
    state_ok = true;
    return 0;
  }

  int notify(char *ptr){
    if( !state_ok ){
      return -1;
    }
    size_t size = strlen(ptr) + 1 /* '\0' */;
    if( size > MAX_FILE_PATH ){
      printf("ckpt_compressor: file name is too big, max is %d\n", MAX_FILE_PATH);
      return -1;
    }
    int ret;
    if( (ret = write(wfd,ptr,size) ) < 0 ){
      perror("Compressor write");
      state_ok = false;
      return ret;
    }
    return 0;
  }

  int finish(){
    // compressor_thread will notice that fd was closed and exit
    if( state_ok ){
      close(wfd);
      int status;
      wait(&status);
    }
    return 0;
  }
};

#endif // COMPRESS_H
