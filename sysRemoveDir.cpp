#include <sys/types.h>
#include <dirent.h>
#include <string>
#include <cstring>
#include <stdio.h>
#include <unistd.h>

void sysRemoveDir(std::string dirName) {//remove directory and all its contents
    DIR *dirp = opendir(dirName.c_str());
    if (dirp != NULL) {//otherwise directory does not exist
        struct dirent *dp = NULL;
        while ((dp = readdir(dirp)) != NULL) {
            std::string file1=dp->d_name;
            if (file1!="." && file1!="..") {
                remove((dirName+file1).c_str());
            };
        };
        (void)closedir(dirp);
        rmdir(dirName.c_str());
    };
};
