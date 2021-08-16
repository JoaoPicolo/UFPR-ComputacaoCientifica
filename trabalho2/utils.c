// ###################################
// Gabriel Marczuk Thá  - GRR20186070 
// João Pedro Picolo    - GRR20182659
// ###################################

#include "utils.h"
#include <sys/time.h>

/*!
    \brief Pega tempo da atual da máquina.

    \return Tempo atual da máquina em milissegundos.
*/
double timestamp() {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}
