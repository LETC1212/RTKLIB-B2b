#include "rtapp.h"
#include "rtklib.h"

PPPGlobal_t PPP_Glo = {0};

int main(int argc, char **argv)
{
    //app_convbin(argc, argv);

    app_rtkrcv(argc, argv);

    return 0;
}