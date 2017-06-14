#ifndef __BLOCK_FILE
#define __BLOCK_FILE

//header files
#include "HRTree.h"

class BlockFile{
private:
    FILE* file_ptr;     // os file pointer
    char* file_name;    // os file name
    int block_len;      // length of a block
    int block_pos;      // block # of fp's position (fp can stay at block boundaries)
    int block_num;      // total # of blocks
    bool new_flag;      // specifies if this is a new file

public:
    BlockFile(const char* filename, int block_len);
    ~BlockFile();

private:
    void put_bytes(char* bytes, int num){fwrite(bytes, num, 1, file_ptr);}
    void get_bytes(char* bytes, int num){fread(bytes, num, 1, file_ptr);}
    void fwrite_num(int num);
    int fread_num();
public:
    void seek_block(int block_num){fseek(file_ptr,(block_num - block_pos)*blocklength,SEEK_CUR);}
    void read_block(Block b, int index);
    void write_blodk(Block b, int index);

    int append_block(Block b);
    bool delete_last_block(int num);
    bool file_new(){return new_flag;}
    int get_block_length(){return block_len;}
    int get_num_of_blocks(){return block_num;}
};

#endif //__BLOCK_FILE
