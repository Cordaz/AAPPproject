
unsigned long djb2_hash(unsigned char *str);

uint64_t MurmurHash64A (const void * key, int len, unsigned int seed);

uint64_t APHash(char* str, unsigned int length);

uint64_t fnvhash(char * string);

uint64_t SDBMHash(char* str, unsigned int length);

uint64_t RSHash(char* str, unsigned int length);

