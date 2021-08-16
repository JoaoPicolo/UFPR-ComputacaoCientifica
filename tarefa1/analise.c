#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "types.h"
#include "operations.h"
#include "manipulation.h"

void saveInput(char firstOperator[], char secondOperator[],
                char operator, operation *operation) {
    // Removes X from string in order to get index
    memmove(firstOperator, firstOperator + 1, strlen(firstOperator));
    memmove(secondOperator, secondOperator + 1, strlen(secondOperator));

    // Save infos converting to long int, default type for strtol function
    // Less 1 to reference intervals array
    operation->firstIndex = strtol(firstOperator, (char **)NULL, 10) - 1;
    operation->secondIndex = strtol(secondOperator, (char **)NULL, 10) - 1;
    operation->op= operator;
}

int main() {
    int m, n;
    interval *intervals;
    operation *operations;

    scanf("%d %d", &m, &n);

    intervals = malloc((m + n) * sizeof(interval));
    operations = malloc(n * sizeof(operation));

    
    Float_t num;
    float min, max;
    char input[10];
    for(int i = 0; i < m; i++) {
        scanf("%s %f", input, &num.f);
        min = getMinEqual(num);
        max = getMaxEqual(num);
        if(!isnan(min) && !isnan(max)) {
            intervals[i] = getInterval(min, max);
        }
        else {
            if(isnan(min)) {
                intervals[i] = getInterval(-1.0 * max, max);
            }
            else {
                intervals[i] = getInterval(-1.0 * min, min);
            }
        }

        intervals[i].isUnitary = isUnitary(intervals[i].minimum, intervals[i].maximum) ? 1 : 0;
    }

    char firstOperator[10], secondOperator[10], operator;
    for(int i = 0; i < n; i++) {
	    scanf("%s %s %s %c %s", input, input, firstOperator, &operator, secondOperator);
	    saveInput(firstOperator, secondOperator, operator, &operations[i]);
    }

    int isValid = 1;
    interval firstInterval, secondInterval;
    for(int i = 0; i < n; i++) {
        operator = operations[i].op;
        firstInterval = intervals[operations[i].firstIndex];
        secondInterval = intervals[operations[i].secondIndex];
        switch(operator) {
            case '+':
                isValid = sumInterval(firstInterval, secondInterval, &intervals[i + m]);
                if(isValid == -1) {
                    return isValid;
                }
                break;
            
            case '-': 
                isValid = subtractInterval(firstInterval, secondInterval, &intervals[i + m]);
                if(isValid == -1) {
                    return isValid;
                }
                break;

            case '*': 
                isValid = multiplyInterval(firstInterval, secondInterval, &intervals[i + m]);
                if(isValid == -1) {
                    return isValid;
                }
                break;

            case '/': 
                isValid = divideInterval(firstInterval, secondInterval, &intervals[i + m]);
                if(isValid == -1) {
                    return isValid;
                }
                break;

            default:
                printf("Undefined operation\n");
                return 0;
        }
    }


    // Print data
    for(int i = 0; i < (m + n); i++) {
	    printf("X%d = [%1.8e, %1.8e]\n", i + 1, intervals[i].minimum, intervals[i].maximum);
    }

    printf("Não unitários:\n");
    for(int i = m; i < m + n; i++) {
        if(!intervals[i].isUnitary) {
            printf("X%d = [%1.8e, %1.8e]\n", i + 1, intervals[i].minimum, intervals[i].maximum);
        }
    }

    free(intervals);
    free(operations);
    return 0;
}
