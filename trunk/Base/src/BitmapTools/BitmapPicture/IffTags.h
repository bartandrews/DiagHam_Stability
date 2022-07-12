

#include "config.h"


#ifdef __LITTLEENDIAN__

// Bitmap Header tag
#define BMHD_TAG 0x424d4844
// Body tag
#define BODY_TAG 0x424f4459
// Color Map tag
#define CMAP_TAG 0x434d4150
// IFF Format File tag
#define FORM_TAG 0x464f524d
// IFF Bitmap File Format tag
#define ILBM_TAG 0x494c424d


#else

// Bitmap Header tag
#define BMHD_TAG 0x44484d42
// Body tag
#define BODY_TAG 0x59444f42
// Color Map tag
#define CMAP_TAG 0x50414d43
// IFF Format File tag
#define FORM_TAG 0x4d524f46
// IFF Bitmap File Format tag
#define ILBM_TAG 0x4d424c49

#endif
