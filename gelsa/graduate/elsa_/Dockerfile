FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget sudo git gcc g++ python3.8 python3.8-dev python3-pip\
    &&sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1\
    &&python -m pip install numpy pandas graphistry scipy

RUN git clone https://github.com/pybind/pybind11.git
RUN git clone https://github.com/galaxyproject/galaxy.git
RUN cd /galaxy/
# RUN git checkout bb8c154899dbfae4f84f1b2a1f089bfd3ec4aedd

# docker build -t galaxy_elsa:1.0 .
# docker run -d --name mycon galaxy_elsa:1.0  tail -f /dev/null
# docker tag galaxy_elsa:1.1 bingwujiayin/galaxy_elsa:1.0

# docker——environment
sudo apt update
sudo apt install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io
docker images

sudo apt-get update && apt-get install -y wget \
&&wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.1-1_all.deb \
&&dpkg -i cuda-keyring_1.1-1_all.deb && apt-get update \
&&apt-get install -y cuda 

PATH="/usr/local/cuda-12.4/bin:${PATH}"
LD_LIBRARY_PATH="/usr/local/cuda-12.4/lib64:${LD_LIBRARY_PATH}"


# docker run --gpus all --name gu ubuntu:jammy-20240227 tail -f /dev/null

