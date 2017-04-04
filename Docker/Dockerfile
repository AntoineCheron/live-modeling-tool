FROM ubuntu
LABEL name="live-modeling-tool"
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y unzip
RUN apt-get install -y curl
RUN curl -sL https://deb.nodesource.com/setup_7.x -o nodesource_setup.sh 
RUN chmod +x nodesource_setup.sh
RUN ./nodesource_setup.sh
RUN apt-get install -y nodejs
RUN apt-get install -y git
RUN apt-get install -y openjdk-8-jdk
RUN apt-cache search jdk
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
RUN export PATH=$PATH:/usr/lib/jvm/java-8-openjdk-amd64/bin
RUN apt install -y maven
RUN git clone https://github.com/AntoineCheron/live-modeling-tool.git
WORKDIR live-modeling-tool/
RUN git fetch
RUN git checkout releases
RUN git pull
RUN ./install.sh
COPY startup.sh /
WORKDIR /
RUN chmod +x /startup.sh
ENTRYPOINT ["./startup.sh"]


