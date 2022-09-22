FROM ubuntu:20.04

RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list \
&& sed -i s@/security.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list \
&& apt-get clean && apt-get update

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ Asia/Shanghai
# ENV LANG zh_CN.UTF-8

RUN apt-get install --no-install-recommends -y build-essential python3.8 python3.8-dev python3-pip octave\
 && ln -s /usr/bin/python3.8 /usr/bin/python \
 && pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple \
 && apt-get clean \
 && rm -rf /tmp/* /var/lib/apt/lists/* /var/tmp/*

RUN pip install numpy scipy wandb tqdm\
 && rm -rf ~/.cache/pip

WORKDIR /home/svd_competition