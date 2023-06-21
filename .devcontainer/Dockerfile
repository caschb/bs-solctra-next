FROM archlinux:base-devel

RUN pacman --noconfirm -Syu && pacman -S --noconfirm clang openmpi cmake ninja git
RUN useradd user
USER user