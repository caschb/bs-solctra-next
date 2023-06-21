FROM archlinux:base-devel

ARG USERNAME=myuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN pacman --noconfirm -Syu && pacman -S --noconfirm clang openmp openmpi cmake ninja git
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME
USER $USERNAME