# Live Modeling Tool

Research engineering internship project at IRISA, France.
*Goal :* to develop a live modeling tool for hydrologist researcher. This tool
let them do two main things :
* Run simulation from entry data through a mathematic model. Then, display the
 output data they select onto charts and go through those charts as they need.
 Moreover, let them modify the entry data and see their impact on the results in
 real time
 * Adjust the variables of mathematic model from real data and the adjusting
 algorithm they input into this system.

## Typical Use Cases
* User wants to produce a Boussinesq Simulation on a hillslope and see the results
on a chart, choosing the kind of chart. Then go through the chart.
* User wants to adjust its mathematic model comparing the simulation's results
and the data she picked up on the real case study.

## Getting Started

### Prerequisites
To run this project, you'll need :

*Docker*
```
https://www.docker.com/community-edition
```

*A Java Virtual Machine*
```
https://java.com/en/download/
```

You will also need maven and unzip
```
sudo apt install maven
sudo apt-get install unzip
```

### Installing
You only need to download the jar contained in this project in ./Jar/.
First : open your terminal and go to the folder where you want to install this project

Then,
For the stable version (latest release) :
```
wget "https://github.com/AntoineCheron/live-modeling-tool/archive/v0.1.1.zip"
unzip master.zip -d ./
cd live-modeling-tool-master
mvn package
```

For the last development version :
```
wget "https://github.com/AntoineCheron/live-modeling-tool/archive/master.zip"
unzip master.zip -d ./
cd live-modeling-tool-master
mvn package
```

When you will try to run the simulation, the server will automatically pull the
latest version of the needed docker image, containing the python environment.
If you want to save time during your first execution, you can download the docker image now.
To do that :
```
docker pull antoinecheronirisa/hs1d
```

### Running
1. Go to the folder where you installed the project
2. Run the following command :
```
java -jar ./target/live-modeling-tool-0.1-jar-with-dependencies.jar
```

## How to go through the web app
After having runned the project as described in the Running section
under Getting Started :

1. Open your favorite web browser and go to http://localhost:8080/
2. Input 4 files, containing entry data :
* geologic.input
* hydrologic.input
* morphologic.input
* key_spatialized_param.param

3. Run simulation & wait until it finishes
4. Add a new chart and select the abscissa, the type of chart and the data
you want to see on the ordinate
5. Click the "Generate chart" button to see the chart
6. Enjoy !

## Software Architecture
### Front-end

### Back-end

## Next features that will be developed
* Calibration of the simulation's mathematic model
* Making simulation faster
* Charts displaying matrices

More to come ...

## Technologies

### Front-end
* [Vue.js](https://vuejs.org/) - A vue templating framework
* [UnderscoreJS](http://underscorejs.org/) - A JS library to make arrays and objects manipulation easier
* [NPM](https://www.npmjs.com/) - Web dependencies management
* [Babel](https://babeljs.io/) - A Javascript compiler
* [Webpack](https://webpack.js.org/) - A depencies bundler to improve performances

### Back-end
Developed with Java and [Jetty](http://www.eclipse.org/jetty/)
* [Docker](https://www.docker.com/)
* [Maven](https://maven.apache.org/) - Dependency management
* [Python](https://www.python.org/) for the mathematic model used for simulations

## Acknowledgments
* [Weizhenye](https://woozy.im/) for the [vue-highcharts component](https://github.com/weizhenye/vue-highcharts)
* [Bootswatch](https://bootswatch.com/solar/) for the Solar CSS Bootstrap theme
* [Benjamin Fox](https://github.com/foxbenjaminfox) for the [vue-async-computed library](https://github.com/foxbenjaminfox/vue-async-computed)
