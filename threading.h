#include <pthread.h>
#include <string.h>
#include <stdint.h>

// Spinlock-based threading implementation

const size_t NUM_THREADS = 8;
const size_t SIGNAL_SIZE = (NUM_THREADS - 1) / sizeof(uint64_t) + 1;

typedef struct {
    pthread_t threads[NUM_THREADS];
    double *thread_accelerations[NUM_THREADS];
    volatile uint64_t signals[SIGNAL_SIZE];
    volatile uint8_t flag_signal;
    volatile bool die_flag;
} thread_handler_t;

const uint8_t START_SIGNAL = 2;
const uint8_t COLLECT_SIGNAL = 1;
const uint8_t DONE_SIGNAL = 0;

static inline uint64_t repeat_byte(uint8_t byte, size_t n) {
    uint64_t out = byte;
    for (size_t i = 1; i < n; i++) {
        out |= ((uint64_t) byte) << (8 * i);
    }
    return out;
}

bool check_signals(thread_handler_t *handler, uint8_t signal) {
    uint64_t check = repeat_byte(signal, sizeof(uint64_t));
    for (size_t i = 0; i < NUM_THREADS/8; i++) {
        if (handler->signals[i] != check) {
            return false;
        }
    }
    const size_t DIFF = SIGNAL_SIZE * sizeof(uint64_t) - NUM_THREADS;
    if (DIFF) {
        return handler->signals[SIGNAL_SIZE - 1] == repeat_byte(signal, DIFF);
    } else {
        return true;
    }
}

void wait_signals(thread_handler_t *handler, uint8_t signal) {
    while (!check_signals(handler, signal))
    ;
}

void wait_flag(thread_handler_t *handler, uint8_t signal) {
    while (!handler->die_flag && handler->flag_signal != signal)
    ;
}

void set_all_signals(thread_handler_t *handler, uint8_t signal) {
    for (size_t i = 0; i < NUM_THREADS/8; i++) {
        handler->signals[i] = repeat_byte(signal, sizeof(uint64_t));
    }
    const size_t DIFF = SIGNAL_SIZE * sizeof(uint64_t) - NUM_THREADS;
    if (DIFF) {
        handler->signals[SIGNAL_SIZE - 1] = repeat_byte(signal, DIFF);
    }
}

void set_signal(thread_handler_t *handler, uint8_t signal, size_t tid) {
    ((volatile uint8_t *) handler->signals)[tid] = signal;
}

void set_flag(thread_handler_t *handler, uint8_t signal) {
    handler->flag_signal = signal;
}