#include <stdio.h>
#include <stdlib.h>
#include "gifer.h"

/*
 * Simple GIF encoding, implementation borrows VERY heavily
 * from lecram's gifenc library: https://github.com/lecram/gifenc
 */

#define HAS_GCT 128

#define lo_byte(n) ((n) & 0xff)
#define hi_byte(n) ((n) >> 8)

unsigned int bitlen(unsigned int x)
{
    unsigned int l = 0;
    while(x > 0)
    {
        x >>= 1;
        l++;
    }

    return l;
}

unsigned int pow2(unsigned int x)
{
    unsigned int y = 1;
    while(x > 0)
    {
        y <<= 1;
        x--;
    }

    return y;
}

void init_table(struct node *n, uint16_t num_colors)
{
    n->code = UINT16_MAX;
    n->num_colors = num_colors;
    n->children = (struct node *)malloc(num_colors * sizeof(struct node));

    int i;
    for(i=0; i<num_colors; i++)
    {
        n->children[i].code = i;
        n->children[i].num_colors = num_colors;
        n->children[i].children = NULL;
    }
}

void init_node(struct node *n)
{
    n->children = (struct node *)malloc(n->num_colors * sizeof(struct node));

    int i;
    for(i=0; i<n->num_colors; i++)
    {
        n->children[i].code = UINT16_MAX;
        n->children[i].num_colors = n->num_colors;
        n->children[i].children = NULL;
    }
}

void free_node(struct node *n)
{
    int i;
    for(i=0; i<n->num_colors; i++)
        if(n->children[i].children != NULL)
            free_node(&(n->children[i]));
    free(n->children);
    n->children = NULL;
    n->num_colors = 0;
    n->code = 0;
}


void update_block(uint8_t *block, uint32_t blocksize, uint16_t code,
                  uint8_t codesize, uint16_t *idx, uint8_t *offset,
                  uint32_t *buffer, FILE *f)
{
    uint32_t mask = (1u << codesize) - 1u;  //make mask for codesize
    *buffer &= ~(mask << *offset);          //0 dest bits in buffer
    *buffer |= ((mask & code) << *offset);  //set buffer bits with code
    *offset += codesize;

    while(*offset >= 8)
    {
        block[*idx] = (uint8_t) *buffer;  //Grabs the least significant byte
        *idx += 1;
        *offset -= 8;
        *buffer >>= 8;

        if(*idx == blocksize)
        {
            block[0] = 255;
            fwrite(block, sizeof(uint8_t), blocksize, f);

            *idx = 1;
        }
    }

    block[*idx] = (uint8_t) *buffer;
}

void flush_block(uint8_t *block, uint32_t blocksize, uint16_t idx,
                 uint8_t offset, FILE *f)
{
    if(idx > 1 || offset > 0)
    {
        block[0] = idx;
        fwrite(block, sizeof(uint8_t), idx+1, f);
    }

    block[0] = 0;
    fwrite(block, sizeof(uint8_t), 1, f);
}

void print_code_bits(uint8_t code, uint8_t codesize)
{
    char s[codesize+1];
    int i;
    uint8_t mask = 1u;
    for(i=codesize-1; i >= 0; i--)
    {
        s[i] = (mask & code) ? '1' : '0';
        mask <<= 1;
    }
    s[codesize] = '\0';

    printf("%s", s);
}

void print_buffer_bits(uint32_t buffer)
{
    print_code_bits(buffer >> 24, 8);
    printf(" ");
    print_code_bits(buffer >> 16, 8);
    printf(" ");
    print_code_bits(buffer >> 8, 8);
    printf(" ");
    print_code_bits(buffer, 8);
}

void writeFrame(uint16_t w, uint16_t h, int *data, int idx_tp,
                uint16_t gct_size_code, size_t gctsize, uint16_t delaytime,
                FILE *f)
{
    // Graphics Control Extension
    
    int use_tp = (idx_tp >= 0 && idx_tp < 256) ? 1 : 0;

    uint8_t gce[8];
    gce[0] = 0x21;  // extension block
    gce[1] = 0xf9;  // graphics extension
    gce[2] = 0x04;  // 4 bytes of data follow
    gce[3] = 0x00 | (use_tp ? 0x01 : 0x00) | 0x08;
    gce[4] = lo_byte(delaytime);  //Delay Time lo
    gce[5] = hi_byte(delaytime);  //Delay Time hi
    gce[6] = use_tp ? (uint8_t) idx_tp : 0;
    gce[7] = 0x00;

    fwrite(gce, sizeof(uint8_t), 8, f);

    // Image Descriptor

    uint8_t id[10];
    id[0] = 0x2c;   //Image separator, always the first byte.
    id[1] = 0x00;   //Left edge
    id[2] = 0x00;   //
    id[3] = 0x00;   //Top edge
    id[4] = 0x00;   //
    id[5] = lo_byte(w);   //Width
    id[6] = hi_byte(w);   //
    id[7] = lo_byte(h);   //Height
    id[8] = hi_byte(h);   //
    id[9] = 0x00;   // Other things

    fwrite(id, sizeof(uint8_t), 10, f);

    // Image Data

    uint8_t min_lzw_code_size = gct_size_code+1;
    fwrite(&min_lzw_code_size, sizeof(uint8_t), 1, f);

    uint8_t block[256];
    block[0] = 0;
    uint8_t offset = 0;
    uint16_t idx = 1;
    uint32_t buffer = 0;

    uint8_t codesize = min_lzw_code_size + 1;
    uint16_t clear_code = gctsize;
    uint16_t eoi_code = gctsize + 1;
    uint16_t next_code = gctsize + 2;

    struct node root;
    init_table(&root, gctsize);
    int npix = ((int) w) * ((int) h);

    struct node *entry = &root;
    /*
    printf("Emit code: ");
    print_code_bits(clear_code, codesize);
    printf("(%d)\n", clear_code);

    printf("Buffer: ");
    print_buffer_bits(buffer);
    printf("\n");
    */

    update_block(block, 256, clear_code, codesize, &idx, &offset, &buffer, f);

    /*
    printf("Buffer: ");
    print_buffer_bits(buffer);
    printf("\n");
    */

    int i;

    for(i=0; i < npix; i++)
    {
        uint8_t pix = data[i];
        //printf("%d | %d %d %d %d\n", i, pix, idx, offset, codesize);

        if(entry->children && entry->children[pix].code < next_code)
            entry = &(entry->children[pix]);
        else
        {
            update_block(block, 256, entry->code, codesize, &idx, &offset,
                         &buffer, f);

            /*
            printf("Buffer: ");
            print_buffer_bits(buffer);
            printf("\n");
            */

            if(entry->children == NULL)
                init_node(entry);
            (entry->children[pix]).code = next_code;

            /*
            printf("Emit code: ");
            print_code_bits(entry->code, codesize);
            */
        
            if(next_code & (1u << codesize))
                codesize++;

            /*
            printf("(%d) Add code: ", entry->code);
            print_code_bits(next_code, codesize);
            printf(" (%d)\n", next_code);
            */

            if(codesize >= 13)
            {
                update_block(block, 256, clear_code, 12,
                             &idx, &offset, &buffer, f);
                free_node(&root);
                init_table(&root, gctsize);
                codesize = min_lzw_code_size + 1;
                next_code = gctsize + 2;
            }
            else
                next_code++;

            entry = &(root.children[pix]);
        }

    }
    update_block(block, 256, entry->code, codesize, &idx, &offset, &buffer, f);
    update_block(block, 256, eoi_code, codesize, &idx, &offset, &buffer, f);

    flush_block(block, 256, idx, offset, f);

    free_node(&root);
}

void makeGIF(int width, int height, int Nframes, int *palette, int Nc,
             int idx_bg, int idx_tp, double framerate, int *data,
             const char *filename)
{
    if(width >= 65536 || height >= 65536)
    {
        printf("Image too big or small for GIF\n");
        return;
    }
    if(width < 0 || height < 0)
    {
        printf("Image dimensions negative\n");
        return;
    }
    if(Nc > 256)
    {
        printf("Too many colors for GIF!\n");
        return;
    }
    if(Nc <= 1)
    {
        printf("Not enough colors for GIF!\n");
        return;
    }

    uint16_t w = (uint16_t) width;
    uint16_t h = (uint16_t) height;
    uint16_t gct_size_code = bitlen(Nc-1)-1;

    size_t gctsize = pow2(gct_size_code+1);

    FILE *f = fopen(filename, "wb");

    // Header
    char hdr[] = "GIF89a";
    fwrite(hdr, sizeof(char), 6, f);

    // Logical Screen Descriptor

    uint8_t lsd[7];
    lsd[0] = lo_byte(w);
    lsd[1] = hi_byte(w);
    lsd[2] = lo_byte(h);
    lsd[3] = hi_byte(h);
    lsd[4] = 0xF0 | gct_size_code;
    lsd[5] = (idx_bg >= 0 && idx_bg < 256) ? (uint8_t) idx_bg : 0;
    lsd[6] = 0;

    fwrite(lsd, sizeof(uint8_t), 7, f);

    //Global Color Table
    int i;

    uint8_t gct[3*gctsize];
    for(i=0; i<3*Nc; i++)
        gct[i] = (uint8_t) palette[i];
    for(i=3*Nc; i<3*gctsize; i++)
        gct[i] = 0;

    fwrite(gct, sizeof(uint8_t), 3*gctsize, f);

    //Write the image frames
    if(Nframes > 1)
    {
        // Application Extension
        uint8_t appext_hdr[3];
        appext_hdr[0] = 0x21;  // Extension
        appext_hdr[1] = 0xFF;  // Application Extension
        appext_hdr[2] = 11u;    // 11 byte application block
        unsigned char appext_blk[11] = "NETSCAPE2.0";
        uint8_t appext_subblk[5];
        appext_subblk[0] = 3u;
        appext_subblk[1] = 1u;
        appext_subblk[2] = 0;
        appext_subblk[3] = 0;
        appext_subblk[4] = 0;
        fwrite(appext_hdr, sizeof(uint8_t), 3, f);
        fwrite(appext_blk, sizeof(unsigned char), 11, f);
        fwrite(appext_subblk, sizeof(uint8_t), 5, f);

        int i;
        size_t npix = width * height;

        int delaytime_int = 100.0 / framerate;
        uint16_t delaytime;
        if(delaytime_int <= (uint16_t) 0xFFFF)
            delaytime = (uint16_t) delaytime_int;
        else
            delaytime = (uint16_t) 0xFFFF;

        for(i=0; i<Nframes; i++)
            writeFrame(w, h, data + npix*i, idx_tp, gct_size_code, gctsize, 
                       delaytime, f);
    }
    else
        writeFrame(w, h, data, idx_tp, gct_size_code, gctsize, 0, f);

    //trailer
    
    uint8_t trailer = 0x3b;
    fwrite(&trailer, sizeof(uint8_t), 1, f);

    fclose(f);
}
