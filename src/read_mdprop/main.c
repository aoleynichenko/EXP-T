#include <complex.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "mdprop.h"
#include "mrconee.h"
#include "unformatted.h"


int main()
{
    mrconee_data_t *mrconee_data = read_mrconee("MRCONEE");
    if (mrconee_data == NULL) {
        printf(" MRCONEE file not found\n");
    }
    else {
        print_mrconee_data(stdout, mrconee_data);
    }

    read_mdprop("MDPROP", mrconee_data);

    return 0;
}



