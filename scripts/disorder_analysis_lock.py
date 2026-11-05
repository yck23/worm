"""Cross-process locking for shared notebook/plot refreshes."""
import errno
import os
from pathlib import Path
import time


class AnalysisLock:
    def __init__(self,path,timeout=540):
        self.path=Path(path);self.timeout=timeout;self.stream=None

    def acquire(self):
        self.path.parent.mkdir(parents=True,exist_ok=True)
        self.stream=self.path.open('a+b')
        if self.path.stat().st_size==0:self.stream.write(b'0');self.stream.flush()
        deadline=time.monotonic()+self.timeout
        while True:
            self.stream.seek(0)
            try:
                if os.name=='nt':
                    import msvcrt
                    msvcrt.locking(self.stream.fileno(),msvcrt.LK_NBLCK,1)
                else:
                    import fcntl
                    fcntl.flock(self.stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
                return self
            except OSError as exc:
                if exc.errno not in (errno.EACCES,errno.EAGAIN,errno.EDEADLK):
                    self.stream.close();self.stream=None;raise
                if time.monotonic()>=deadline:
                    self.stream.close();self.stream=None
                    raise TimeoutError('Another process is refreshing disorder figures') from exc
                time.sleep(.1)

    def release(self):
        if self.stream is None:return
        try:
            self.stream.seek(0)
            if os.name=='nt':
                import msvcrt
                msvcrt.locking(self.stream.fileno(),msvcrt.LK_UNLCK,1)
            else:
                import fcntl
                fcntl.flock(self.stream.fileno(),fcntl.LOCK_UN)
        finally:self.stream.close();self.stream=None
