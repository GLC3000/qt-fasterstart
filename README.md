# qt-fasterstart
Instantly fast-start MP4 with filesystem block cloning

---

This is a fork of [ffmpeg's qt-faststart](https://github.com/FFmpeg/FFmpeg/blob/master/tools/qt-faststart.c).
It rearranges an MP4 file such that the moov atom is in front of the data.
Instead of copying the data, it uses filesystem block cloning API provided by [Linux](https://www.man7.org/linux/man-pages/man2/ioctl_ficlonerange.2.html) or [Windows](https://learn.microsoft.com/en-us/windows/win32/fileio/block-cloning) to save time as well as disk.

At the time of writing, it supports [XFS](https://blogs.oracle.com/linux/post/xfs-data-block-sharing-reflink) and Btrfs on Linux. It also supports [ReFS](https://learn.microsoft.com/en-us/windows-server/storage/refs/refs-overview) and [Btrfs on Windows](https://github.com/maharmstone/btrfs).

To satisfy the filesystem block alignment requirement, it inserts a small (usually less than 4096 bytes) padding area between the moov atom and the rest of data.
