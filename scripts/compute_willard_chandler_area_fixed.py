#!/usr/bin/env python3
"""Example showing proper timing for both CPU and GPU backends."""

import time
from typing import Optional

try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None


class Timer:
    """Backend-aware timer that uses appropriate timing method."""

    def __init__(self, backend: str):
        self.backend = backend.lower()
        self.use_cuda_events = (self.backend == 'cupy' and CUPY_AVAILABLE)

        if self.use_cuda_events:
            self.start_event = cp.cuda.Event()
            self.end_event = cp.cuda.Event()

        self.count = 0
        self.total_time = 0.0

    def start(self):
        """Start timing."""
        if self.use_cuda_events:
            self.start_event.record()
        else:
            self.start_time = time.perf_counter()

    def stop(self):
        """Stop timing and accumulate."""
        if self.use_cuda_events:
            self.end_event.record()
            self.end_event.synchronize()
            elapsed_ms = cp.cuda.get_elapsed_time(self.start_event, self.end_event)
        else:
            elapsed_ms = (time.perf_counter() - self.start_time) * 1000.0

        self.count += 1
        self.total_time += elapsed_ms

    def get_average_ms(self) -> Optional[float]:
        """Get average time in milliseconds."""
        if self.count == 0:
            return None
        return self.total_time / self.count

    def get_total_ms(self) -> float:
        """Get total time in milliseconds."""
        return self.total_time


# Example usage in the main loop:
def example_usage(backend='cpu'):
    timer = Timer(backend)

    # In your frame loop:
    for frame in range(10):
        timer.start()

        # ... do WillardChandler computation ...

        timer.stop()

    avg_time = timer.get_average_ms()
    print(f"Average time per frame: {avg_time:.2f} ms")
    print(f"Total time: {timer.get_total_ms():.2f} ms")


if __name__ == "__main__":
    print("CPU timing:")
    example_usage(backend='cpu')

    if CUPY_AVAILABLE:
        print("\nGPU timing:")
        example_usage(backend='cupy')
