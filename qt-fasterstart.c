/*
 * qt-fasterstart.c, v0.2
 * adapted from ffmpeg version by Mike Melanson
 * This file is placed in the public domain. Use the program however you
 * see fit.
 *
 * This utility rearranges a Quicktime file such that the moov atom
 * is in front of the data, thus facilitating network streaming.
 *
 * Compile the program with:
 *  gcc -O2 -o qt-fasterstart qt-fasterstart.c
 *
 * Invoke the program with:
 *  qt-fasterstart <infile.mov> <outfile.mov>
 *
 * Notes: Quicktime files can come in many configurations of top-level
 * atoms. This utility stipulates that the very last atom in the file needs
 * to be a moov atom. When given such a file, this utility will rearrange
 * the top-level atoms by shifting the moov atom from the back of the file
 * to the front, and patch the chunk offsets along the way. This utility
 * presently only operates on uncompressed moov atoms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <limits.h>

#ifdef _WIN32
#undef fseeko
#define fseeko(x, y, z) _fseeki64(x, y, z)
#undef ftello
#define ftello(x)       _ftelli64(x)
#define MIN(a,b) ((a) > (b) ? (b) : (a))

#include <windows.h>
#define FSCTL_DUPLICATE_EXTENTS_TO_FILE 0x00098344
typedef struct _DUPLICATE_EXTENTS_DATA {
    HANDLE  FileHandle;
    int64_t SourceFileOffset;
    int64_t TargetFileOffset;
    int64_t ByteCount;
} DUPLICATE_EXTENTS_DATA;

HANDLE _get_osfhandle(int fd);
__int64 _lseeki64(int fd, __int64 offset, int origin);

int clone(int src, int dst, uint64_t src_offset, uint64_t dst_offset, uint64_t len)
{
    DWORD _;
    HANDLE hIn = _get_osfhandle(src);
    HANDLE hOut = _get_osfhandle(dst);
    DUPLICATE_EXTENTS_DATA data;
    uint64_t bytes;

    while (len) {
        bytes = MIN(len, 1ULL<<30); // 1GB at a time
        data.FileHandle = hIn;
        data.SourceFileOffset = src_offset;
        data.TargetFileOffset = dst_offset;
        data.ByteCount = bytes;
        if(!DeviceIoControl(hOut, FSCTL_DUPLICATE_EXTENTS_TO_FILE, &data, sizeof data, NULL, 0, &_, NULL)) return 0;
        src_offset += bytes;
        dst_offset += bytes;
        len -= bytes;
    }
    return 1;
}

int setsize(int fd, int64_t len)
{
    _lseeki64(fd, len, SEEK_SET);
    return SetEndOfFile(_get_osfhandle(fd));
}
#elif __linux__
#include <unistd.h>
#include <linux/fs.h>
#include <sys/ioctl.h>
#include <sys/statvfs.h>

int clone(int src, int dst, uint64_t src_offset, uint64_t dst_offset, uint64_t len)
{
    struct file_clone_range clone_range;
    clone_range.src_fd = src;
    clone_range.src_offset = src_offset;
    clone_range.src_length = len;
    clone_range.dest_offset = dst_offset;
    return ioctl(dst, FICLONERANGE, &clone_range) != -1;
}

int setsize(int fd, int64_t len)
{
    return ftruncate(fd, len) == 0;
}
#else
#error "Unsupported OS"
#endif

#define ALIGN(x, a)     ALIGN_MASK(x, (a) - 1)
#define ALIGN_MASK(x, mask) (((x) + (mask)) & ~(mask))
#define ALIGN_DOWN(x, a)    ALIGN((x) - ((a) - 1), (a))

#define BE_32(x) (((uint32_t)(((uint8_t*)(x))[0]) << 24) |  \
                             (((uint8_t*)(x))[1]  << 16) |  \
                             (((uint8_t*)(x))[2]  <<  8) |  \
                              ((uint8_t*)(x))[3])

#define BE_64(x) (((uint64_t)(((uint8_t*)(x))[0]) << 56) |  \
                  ((uint64_t)(((uint8_t*)(x))[1]) << 48) |  \
                  ((uint64_t)(((uint8_t*)(x))[2]) << 40) |  \
                  ((uint64_t)(((uint8_t*)(x))[3]) << 32) |  \
                  ((uint64_t)(((uint8_t*)(x))[4]) << 24) |  \
                  ((uint64_t)(((uint8_t*)(x))[5]) << 16) |  \
                  ((uint64_t)(((uint8_t*)(x))[6]) <<  8) |  \
                  ((uint64_t)( (uint8_t*)(x))[7]))

#define AV_WB32(p, val)    {                    \
    ((uint8_t*)(p))[0] = ((val) >> 24) & 0xff;  \
    ((uint8_t*)(p))[1] = ((val) >> 16) & 0xff;  \
    ((uint8_t*)(p))[2] = ((val) >> 8) & 0xff;   \
    ((uint8_t*)(p))[3] = (val) & 0xff;          \
    }

#define AV_WB64(p, val)    {                    \
    AV_WB32(p, (val) >> 32)                     \
    AV_WB32(p + 4, val)                         \
    }

#define BE_FOURCC(ch0, ch1, ch2, ch3)           \
    ( (uint32_t)(unsigned char)(ch3)        |   \
     ((uint32_t)(unsigned char)(ch2) <<  8) |   \
     ((uint32_t)(unsigned char)(ch1) << 16) |   \
     ((uint32_t)(unsigned char)(ch0) << 24) )

#define QT_ATOM BE_FOURCC
/* top level atoms */
#define FREE_ATOM QT_ATOM('f', 'r', 'e', 'e')
#define JUNK_ATOM QT_ATOM('j', 'u', 'n', 'k')
#define MDAT_ATOM QT_ATOM('m', 'd', 'a', 't')
#define MOOV_ATOM QT_ATOM('m', 'o', 'o', 'v')
#define PNOT_ATOM QT_ATOM('p', 'n', 'o', 't')
#define SKIP_ATOM QT_ATOM('s', 'k', 'i', 'p')
#define WIDE_ATOM QT_ATOM('w', 'i', 'd', 'e')
#define PICT_ATOM QT_ATOM('P', 'I', 'C', 'T')
#define FTYP_ATOM QT_ATOM('f', 't', 'y', 'p')
#define UUID_ATOM QT_ATOM('u', 'u', 'i', 'd')

#define CMOV_ATOM QT_ATOM('c', 'm', 'o', 'v')
#define TRAK_ATOM QT_ATOM('t', 'r', 'a', 'k')
#define MDIA_ATOM QT_ATOM('m', 'd', 'i', 'a')
#define MINF_ATOM QT_ATOM('m', 'i', 'n', 'f')
#define STBL_ATOM QT_ATOM('s', 't', 'b', 'l')
#define STCO_ATOM QT_ATOM('s', 't', 'c', 'o')
#define CO64_ATOM QT_ATOM('c', 'o', '6', '4')

#define ATOM_PREAMBLE_SIZE    8
#define COPY_BUFFER_SIZE   33554432
#define MAX_FTYP_ATOM_SIZE 1048576

typedef struct {
    uint32_t type;
    uint32_t header_size;
    uint64_t size;
    unsigned char *data;
} atom_t;

typedef struct {
    uint64_t moov_atom_size;
    uint64_t stco_offset_count;
    uint64_t stco_data_size;
    int stco_overflow;
    uint32_t depth;
} update_chunk_offsets_context_t;

typedef int (*parse_atoms_callback_t)(void *context, atom_t *atom);

static int parse_atoms(
    unsigned char *buf,
    uint64_t size,
    parse_atoms_callback_t callback,
    void *context)
{
    unsigned char *pos = buf;
    unsigned char *end = pos + size;
    atom_t atom;
    int ret;

    while (end - pos >= ATOM_PREAMBLE_SIZE) {
        atom.size = BE_32(pos);
        atom.type = BE_32(pos + 4);
        pos += ATOM_PREAMBLE_SIZE;
        atom.header_size = ATOM_PREAMBLE_SIZE;

        switch (atom.size) {
        case 1:
            if (end - pos < 8) {
                fprintf(stderr, "not enough room for 64 bit atom size\n");
                return -1;
            }

            atom.size = BE_64(pos);
            pos += 8;
            atom.header_size = ATOM_PREAMBLE_SIZE + 8;
            break;

        case 0:
            atom.size = ATOM_PREAMBLE_SIZE + end - pos;
            break;
        }

        if (atom.size < atom.header_size) {
            fprintf(stderr, "atom size %"PRIu64" too small\n", atom.size);
            return -1;
        }

        atom.size -= atom.header_size;

        if (atom.size > end - pos) {
            fprintf(stderr, "atom size %"PRIu64" too big\n", atom.size);
            return -1;
        }

        atom.data = pos;
        ret = callback(context, &atom);
        if (ret < 0) {
            return ret;
        }

        pos += atom.size;
    }

    return 0;
}

static int update_stco_offsets(update_chunk_offsets_context_t *context, atom_t *atom)
{
    uint32_t current_offset;
    uint32_t offset_count;
    unsigned char *pos;
    unsigned char *end;

    printf(" patching stco atom...\n");
    if (atom->size < 8) {
        fprintf(stderr, "stco atom size %"PRIu64" too small\n", atom->size);
        return -1;
    }

    offset_count = BE_32(atom->data + 4);
    if (offset_count > (atom->size - 8) / 4) {
        fprintf(stderr, "stco offset count %"PRIu32" too big\n", offset_count);
        return -1;
    }

    context->stco_offset_count += offset_count;
    context->stco_data_size += atom->size - 8;

    for (pos = atom->data + 8, end = pos + offset_count * 4;
        pos < end;
        pos += 4) {
        current_offset = BE_32(pos);
        if (current_offset > UINT_MAX - context->moov_atom_size) {
            context->stco_overflow = 1;
        }
        current_offset += context->moov_atom_size;
        AV_WB32(pos, current_offset);
    }

    return 0;
}

static int update_co64_offsets(update_chunk_offsets_context_t *context, atom_t *atom)
{
    uint64_t current_offset;
    uint32_t offset_count;
    unsigned char *pos;
    unsigned char *end;

    printf(" patching co64 atom...\n");
    if (atom->size < 8) {
        fprintf(stderr, "co64 atom size %"PRIu64" too small\n", atom->size);
        return -1;
    }

    offset_count = BE_32(atom->data + 4);
    if (offset_count > (atom->size - 8) / 8) {
        fprintf(stderr, "co64 offset count %"PRIu32" too big\n", offset_count);
        return -1;
    }

    for (pos = atom->data + 8, end = pos + offset_count * 8;
        pos < end;
        pos += 8) {
        current_offset = BE_64(pos);
        current_offset += context->moov_atom_size;
        AV_WB64(pos, current_offset);
    }

    return 0;
}

static int update_chunk_offsets_callback(void *ctx, atom_t *atom)
{
    update_chunk_offsets_context_t *context = ctx;
    int ret;

    switch (atom->type) {
    case STCO_ATOM:
        return update_stco_offsets(context, atom);

    case CO64_ATOM:
        return update_co64_offsets(context, atom);

    case MOOV_ATOM:
    case TRAK_ATOM:
    case MDIA_ATOM:
    case MINF_ATOM:
    case STBL_ATOM:
        context->depth++;
        if (context->depth > 10) {
            fprintf(stderr, "atoms too deeply nested\n");
            return -1;
        }

        ret = parse_atoms(
            atom->data,
            atom->size,
            update_chunk_offsets_callback,
            context);
        context->depth--;
        return ret;
    }

    return 0;
}

static int update_moov_atom(
    unsigned char **moov_atom,
    uint64_t *moov_atom_size)
{
    update_chunk_offsets_context_t update_context = { 0 };

    update_context.moov_atom_size = *moov_atom_size;

    if (parse_atoms(
        *moov_atom,
        *moov_atom_size,
        update_chunk_offsets_callback,
        &update_context) < 0) {
        return -1;
    }

    if (!update_context.stco_overflow) {
        return 0;
    }

    fprintf(stderr, "could not upgrade stco atoms to co64...\n");
    return -1;
}

int main(int argc, char *argv[])
{
    FILE *infile  = NULL;
    FILE *outfile = NULL;
    unsigned char atom_bytes[ATOM_PREAMBLE_SIZE];
    uint32_t atom_type   = 0;
    uint64_t atom_size   = 0;
    uint64_t atom_offset = 0;
    int64_t last_offset;
    unsigned char *moov_atom = NULL;
    unsigned char *ftyp_atom = NULL;
    uint64_t moov_atom_size;
    uint64_t free_atom_size;
    uint64_t ftyp_atom_size = 0;
    int64_t start_offset = 0;
    unsigned char *copy_buffer = NULL;
    int bytes_to_copy;
    uint64_t free_size = 0;
    uint64_t moov_size = 0;
    int infd, outfd;
    int64_t current_offset;
    uint64_t clone_length;
    uint64_t fs_block_size = 4096;

    if (argc != 3) {
        printf("Usage: qt-faststart <infile.mov> <outfile.mov>\n"
               "Note: alternatively you can use -movflags +faststart in ffmpeg\n");
        return 0;
    }

    if (!strcmp(argv[1], argv[2])) {
        fprintf(stderr, "input and output files need to be different\n");
        return 1;
    }

    infile = fopen(argv[1], "rb");
    if (!infile) {
        perror(argv[1]);
        goto error_out;
    }
    infd = fileno(infile);

#ifdef __linux__
    struct statvfs vfs;
    fstatvfs(infd, &vfs);
    fs_block_size = vfs.f_bsize;
#endif

    /* traverse through the atoms in the file to make sure that 'moov' is
     * at the end */
    while (!feof(infile)) {
        if (fread(atom_bytes, ATOM_PREAMBLE_SIZE, 1, infile) != 1) {
            break;
        }
        atom_size = BE_32(&atom_bytes[0]);
        atom_type = BE_32(&atom_bytes[4]);

        /* keep ftyp atom */
        if (atom_type == FTYP_ATOM) {
            if (atom_size > MAX_FTYP_ATOM_SIZE) {
                fprintf(stderr, "ftyp atom size %"PRIu64" too big\n",
                       atom_size);
                goto error_out;
            }
            ftyp_atom_size = atom_size;
            free(ftyp_atom);
            ftyp_atom = malloc(ftyp_atom_size);
            if (!ftyp_atom) {
                fprintf(stderr, "could not allocate %"PRIu64" bytes for ftyp atom\n",
                       atom_size);
                goto error_out;
            }
            if (fseeko(infile, -ATOM_PREAMBLE_SIZE, SEEK_CUR) ||
                fread(ftyp_atom, atom_size, 1, infile) != 1 ||
                (start_offset = ftello(infile)) < 0) {
                perror(argv[1]);
                goto error_out;
            }
        } else {
            int ret;
            /* 64-bit special case */
            if (atom_size == 1) {
                if (fread(atom_bytes, ATOM_PREAMBLE_SIZE, 1, infile) != 1) {
                    break;
                }
                atom_size = BE_64(&atom_bytes[0]);
                ret = fseeko(infile, atom_size - ATOM_PREAMBLE_SIZE * 2, SEEK_CUR);
            } else {
                ret = fseeko(infile, atom_size - ATOM_PREAMBLE_SIZE, SEEK_CUR);
            }
            if (ret) {
                perror(argv[1]);
                goto error_out;
            }
        }
        printf("%c%c%c%c %10"PRIu64" %"PRIu64"\n",
               (atom_type >> 24) & 255,
               (atom_type >> 16) & 255,
               (atom_type >>  8) & 255,
               (atom_type >>  0) & 255,
               atom_offset,
               atom_size);
        if ((atom_type != FREE_ATOM) &&
            (atom_type != JUNK_ATOM) &&
            (atom_type != MDAT_ATOM) &&
            (atom_type != MOOV_ATOM) &&
            (atom_type != PNOT_ATOM) &&
            (atom_type != SKIP_ATOM) &&
            (atom_type != WIDE_ATOM) &&
            (atom_type != PICT_ATOM) &&
            (atom_type != UUID_ATOM) &&
            (atom_type != FTYP_ATOM)) {
            fprintf(stderr, "encountered non-QT top-level atom (is this a QuickTime file?)\n");
            break;
        }
        atom_offset += atom_size;

        /* The atom header is 8 (or 16 bytes), if the atom size (which
         * includes these 8 or 16 bytes) is less than that, we won't be
         * able to continue scanning sensibly after this atom, so break. */
        if (atom_size < 8)
            break;

        if (atom_type == MOOV_ATOM)
            moov_size = atom_size;

        if (moov_size && atom_type == FREE_ATOM) {
            free_size += atom_size;
            atom_type = MOOV_ATOM;
            atom_size = moov_size;
        }
    }

    if (atom_type != MOOV_ATOM) {
        printf("last atom in file was not a moov atom\n");
        free(ftyp_atom);
        fclose(infile);
        return 0;
    }

    if (atom_size < 16) {
        fprintf(stderr, "bad moov atom size\n");
        goto error_out;
    }

    /* moov atom was, in fact, the last atom in the chunk; load the whole
     * moov atom */
    if (fseeko(infile, -(atom_size + free_size), SEEK_END)) {
        perror(argv[1]);
        goto error_out;
    }
    last_offset    = ftello(infile);
    if (last_offset < 0) {
        perror(argv[1]);
        goto error_out;
    }

    moov_atom_size = atom_size;
    free_atom_size = fs_block_size - (ftyp_atom_size + moov_atom_size) % fs_block_size;
    if (free_atom_size < 8)
        free_atom_size += fs_block_size;
    free_atom_size += start_offset;
    moov_atom_size += free_atom_size;
    moov_atom      = calloc(1, moov_atom_size - start_offset);
    if (!moov_atom) {
        fprintf(stderr, "could not allocate %"PRIu64" bytes for moov atom\n", atom_size);
        goto error_out;
    }
    if (fread(moov_atom, atom_size, 1, infile) != 1) {
        perror(argv[1]);
        goto error_out;
    }
    AV_WB32(moov_atom + atom_size, free_atom_size);
    AV_WB32(moov_atom + atom_size + 4, FREE_ATOM);

    /* this utility does not support compressed atoms yet, so disqualify
     * files with compressed QT atoms */
    if (BE_32(&moov_atom[12]) == CMOV_ATOM) {
        fprintf(stderr, "this utility does not support compressed moov atoms yet\n");
        goto error_out;
    }

    /* close; will be re-opened later */
    fclose(infile);
    infile = NULL;

    if (update_moov_atom(&moov_atom, &moov_atom_size) < 0) {
        goto error_out;
    }

    /* re-open the input file and open the output file */
    infile = fopen(argv[1], "rb");
    if (!infile) {
        perror(argv[1]);
        goto error_out;
    }
    infd = fileno(infile);

    if (start_offset > 0) { /* seek after ftyp atom */
        if (fseeko(infile, start_offset, SEEK_SET)) {
            perror(argv[1]);
            goto error_out;
        }
    }

    outfile = fopen(argv[2], "wb");
    if (!outfile) {
        perror(argv[2]);
        goto error_out;
    }
    outfd = fileno(outfile);

    /* dump the same ftyp atom */
    if (ftyp_atom_size > 0) {
        printf(" writing ftyp atom...\n");
        if (fwrite(ftyp_atom, ftyp_atom_size, 1, outfile) != 1) {
            perror(argv[2]);
            goto error_out;
        }
    }

    /* dump the new moov atom */
    printf(" writing moov atom...\n");
    if (fwrite(moov_atom, moov_atom_size - start_offset, 1, outfile) != 1) {
        perror(argv[2]);
        goto error_out;
    }

    fflush(infile);
    fflush(outfile);
    current_offset = ftello(outfile);

    clone_length = ALIGN_DOWN(last_offset, fs_block_size);
    printf(" cloning %"PRIu64" bytes...\n", clone_length);

    if (!setsize(outfd, current_offset + last_offset)) {
        fprintf(stderr, "could not set output file size\n");
        goto error_out;
    }

    if (!clone(infd, outfd, 0, current_offset, clone_length)) {
        fprintf(stderr, "could not clone file\n");
#ifdef __linux__
        perror("ioctl_ficlonerange");
#endif
        goto error_out;
    }

    bytes_to_copy = last_offset - clone_length;
    if (bytes_to_copy) {
        fseeko(infile, clone_length, SEEK_SET);
        fseeko(outfile, -bytes_to_copy, SEEK_END);
        copy_buffer = malloc(bytes_to_copy);
        if (!copy_buffer) {
            fprintf(stderr, "could not allocate %d bytes for copy_buffer\n", bytes_to_copy);
            goto error_out;
        }
        printf(" copying rest of file...\n");

        if (fread(copy_buffer, bytes_to_copy, 1, infile) != 1) {
            perror(argv[1]);
            goto error_out;
        }
        if (fwrite(copy_buffer, bytes_to_copy, 1, outfile) != 1) {
            perror(argv[2]);
            goto error_out;
        }
    }

    fclose(infile);
    fclose(outfile);
    free(moov_atom);
    free(ftyp_atom);
    free(copy_buffer);

    return 0;

error_out:
    if (infile)
        fclose(infile);
    if (outfile)
        fclose(outfile);
    free(moov_atom);
    free(ftyp_atom);
    free(copy_buffer);
    return 1;
}
