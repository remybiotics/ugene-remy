/*
 * Bundled htscodecs 1.6.0 declares NEON rANS entry points on aarch64 but does
 * not ship rANS_static32x16pr_neon.c. Map those symbols to the scalar 32x16
 * implementations so Apple Silicon can link.
 */

unsigned char *rans_compress_O0_32x16(unsigned char *in, unsigned int in_size,
                                      unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_32x16(unsigned char *in, unsigned int in_size,
                                        unsigned char *out, unsigned int out_sz);
unsigned char *rans_compress_O1_32x16(unsigned char *in, unsigned int in_size,
                                      unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O1_32x16(unsigned char *in, unsigned int in_size,
                                        unsigned char *out, unsigned int out_sz);

unsigned char *rans_compress_O0_32x16_neon(unsigned char *in, unsigned int in_size,
                                           unsigned char *out, unsigned int *out_size) {
    return rans_compress_O0_32x16(in, in_size, out, out_size);
}

unsigned char *rans_uncompress_O0_32x16_neon(unsigned char *in, unsigned int in_size,
                                             unsigned char *out, unsigned int out_sz) {
    return rans_uncompress_O0_32x16(in, in_size, out, out_sz);
}

unsigned char *rans_compress_O1_32x16_neon(unsigned char *in, unsigned int in_size,
                                           unsigned char *out, unsigned int *out_size) {
    return rans_compress_O1_32x16(in, in_size, out, out_size);
}

unsigned char *rans_uncompress_O1_32x16_neon(unsigned char *in, unsigned int in_size,
                                             unsigned char *out, unsigned int out_sz) {
    return rans_uncompress_O1_32x16(in, in_size, out, out_sz);
}
