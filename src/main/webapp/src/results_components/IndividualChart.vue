<template>
  <div class="generatedChart">
    <highcharts :options="options" v-if="options !== {}"></highcharts>
    <div v-if="matrice">
      <input id="matrice-slider-chart-{{chart.id}}" type="range" class="chart-slider"
        :max="sliderMax"
        @input="handleSliderChange($event.target.value)"
        v-model="sliderValue"
      /> <!-- END OF THE INPUT -->
      <p>
        {{chart.abscissa}} : {{sliderValue}} - Min : 0 - Max : {{sliderMax}}
      </p>
      <button class="btn btn-default" @click="addOneToAbscissa">{{chart.abscissa}} + 1</button>
      <button class="btn btn-default" @click="removeOneToAbscissa">{{chart.abscissa}} - 1</button>
      <button class="btn btn-info" @click="reverseMatrice">
        Reverse matrice</button>
    </div>
  </div>
</template>

<script>
import WS from "../WebServices/WebServices.js"
import Colors from '../utils/colors.js'
import Vue from "vue"
import VueHighcharts from "vue-highcharts"
import AsyncComputed from "vue-async-computed"
import _ from "underscore"

Vue.use(VueHighcharts);
Vue.use(AsyncComputed);

export default {
  name: 'chart',
  props: {
    'chart': {
      type: Object,
      required: true
    }
  },
  data () {
    return {
      options: {},
      matrice: false,
      matriceObject: null,
      reversedMatrice: false,
      sliderMax: 0,
      nbOfMatrices: 0,
      sliderValue: 0,
      changeAbscissaAlert: false,
    }
  },
  asyncComputed: {
    chartWithOptions: function() {
      const parent = this;
      return new Promise(function(resolve, reject) {
        const chart = parent.chart;
        let oneDimensionVariablesFound = false;
        let twoDimensionsVariablesFound = false;
        parent.downloadsDataAndInputThemIntoChart(chart).then(function(newChart) {

          // Makes sure the this.chart object is consistent
          parent.chart = newChart;

          /* Scan the results to see whether there are matrices and/or only
          one dimension variables */
          oneDimensionVariablesFound = parent.containsOneDimensionVariables.apply(this, [chart.selectedResultsValues]);
          twoDimensionsVariablesFound = parent.containsTwoDimensionVariables.apply(this, [chart.selectedResultsValues]);

          /* Prepare the chart object for the generateOptionObjectFromChart function
          depending on the kind of variables */
          if (oneDimensionVariablesFound && twoDimensionsVariablesFound) {
            alert(parent.generateInfoMsgForMultipleNbDimensionsChoice(chart));
            parent.$emit('error');
            reject('Can\'t display one dimension variables and matrices on the same chart');
          } else if (oneDimensionVariablesFound && !twoDimensionsVariablesFound) {
            /* Only one dimension variables to display on the chart. No specific
            computation needed. */
            parent.matrice = false;
          } else if (parent.nbOfMatrices > 1) {
            // Can't display more than one matrice
            alert('Can\'t display more than one matrice per chart');
            parent.$emit('error');
            reject('Can\'t display more than one matrice per chart');
          } else {
            // In this case, the user only wants to display a matrice.
            parent.matrice = true;
            parent.matriceObject = chart.selectedResultsValues[0];

            /* It is possible that the user changed ordinates and/or abscissa
            after having displayed the matrice. In this case, she could have
            reversed the matrice. The following if statement reverse the
            matrice again if the user did it previously. */
            if(parent.reversedMatrice){
              parent.reverseMatrice();
            }

            // The abscissa doesn't change

            // First : place the proper data as the slider
            parent.sliderMax = parent.matriceObject.length;
            // Second : uses matriceObject[0] as the ordinate for the chart
            chart.selectedResultsValues = [parent.matriceObject[0]];
            if(typeof sliderValue !== 'undefined'){
              this.sliderValue = 0;
            }
          }

          // Call this.generateOptionObjectFromChart for the chart to display
          const options = parent.generateOptionObjectFromChart(chart);

          // Finish :D, return the options object for highcharts
          parent.options = options;
          resolve(options);

          // End of parent.downloadsDataAndInputThemIntoChart()
        }).catch(err => {reject(err)});
      });
    }
  },
  methods: {
    handleSliderChange: function(value) {
      // Just need to change the value in the serie of the chart
      const newSerie = this.options.series[0];
      newSerie.data = this.matriceObject[value].map(Number);
      this.options.series[0] = newSerie;
    },
    reverseMatrice: function () {
      if(!this.changeAbscissaAlert) {
        alert('Don\'t forget to change the abscissa');
      }
      // Reverse columns and lines
      const reversedMatrice = [];
      // First : need to initialize a new array containing empty arrays
      for(let i=0; i<this.matriceObject[0].length; i++){
        reversedMatrice[i] = [];
      }
      // Second : acutally reverse the matrice
      for(let i=0; i<this.matriceObject.length; i++){
        for(let j=0; j<this.matriceObject[i].length; j++) {
          reversedMatrice[j][i] = this.matriceObject[i][j];
        }
      }
      // Store the newly calculated matrice and adapt slider
      this.reversedMatrice = !this.reversedMatrice;
      this.matriceObject = reversedMatrice;
      this.sliderValue = 0;
      this.sliderMax = this.matriceObject.length;
      this.handleSliderChange(this.sliderValue);
      this.changeAbscissaAlert = !this.changeAbscissaAlert;
    },
    addOneToAbscissa: function() {
      this.sliderValue ++;
      this.handleSliderChange(this.sliderValue);
    },
    removeOneToAbscissa: function() {
      this.sliderValue --;
      this.handleSliderChange(this.sliderValue);
    },
    generateOptionObjectFromChart: function(chart) {
      // Prepare the xAxis
      const xAxis = chart.abscissaValue.map(Number);

      // Prepare the yAxis and title that will be input into the options object
      // The results (ordinates) are all the data array minus the first element
      const yAxisSeries = [];
      let maxSerieSize = 0;
      let yAxisTitle = "";
      const colors = new Colors();

      for (let i = 0; i < chart.selectedResults.length ; i++) {
        const serieLength = chart.selectedResultsValues[i].length;
        if( serieLength > maxSerieSize){
          maxSerieSize = serieLength;
        }
        // create an yAxis serie correctly formatted for each data
        const tempYAxisSerie = {
          name: chart.selectedResults[i],
          data: chart.selectedResultsValues[i].map(Number),
          color: colors.nextColor(),
          animation: false
        };
        // Add it into the yAxisSeries array
        yAxisSeries.push(tempYAxisSerie);
        // Continue to compute the title
        yAxisTitle += `${chart.selectedResults[i]}, `;
      }
      // Make title looking nicer, removing the coma and space at the end
      yAxisTitle = yAxisTitle.slice(0, yAxisTitle.length -2);

      /* Important step : HERE WE HAVE ALL THE ELEMENTS TO BUILD
      THE OPTIONS OBJECT, LET'S DO THIS*/
      const options = {
        chart: {
          type: chart.type,
          backgroundColor: '#073642'
        },
        title: {
          text: `Chart ${chart.id}`,
          style: {
            color: 'white'
          }
        },
        xAxis: {
          categories: xAxis,
          labels: {
            style: {
              color: 'white'
            }
          }
        },
        yAxis: {
          title: {
            text: yAxisTitle,
            style: {
              color: 'white'
            }
          },
          plotLines: [{
            value: 0,
            width: 1,
            color: '#808080'
          }],
          labels: {
            style: {
              color: 'white'
            }
          }
        },
        legend: {
          layout: 'vertical',
          align : 'right',
          verticalAlign: 'middle',
          borderWidth: 0,
          itemStyle: {
            color: 'white'
          }
        },
        plotOptions: {
          series: {
            turboThreshold: maxSerieSize
          }
        },
        series: yAxisSeries
      };

      /* Finished, now we resolve the promise created on the first line
      of this function */
      return options;
    },
    /* Download the abscissa and the results that the chart need and input them
    all into the chart object */
    downloadsDataAndInputThemIntoChart: function(chart) {
      return new Promise (function(resolve, reject) {
        // First : need to create a Promise Array that will contains all the downloads
        const Downloads = [];

        /* Second : launch every download by pushing a promise into the Downloads array
        for each element to download */
        // 1 : download the abscissa
        Downloads.push(WS.getResult(chart.abscissa));
        // 2 : download all ordinates
        for(let i=0; i < chart.selectedResults.length; i++){
          Downloads.push(WS.getResult(chart.selectedResults[i]));
        }

        /* Third : resolve all the promises, which is the moment when all data
        have been downloaded. Thus, we can check whether there are multiple dimensions
        variables to display or not */
        Promise.all(Downloads).then(function(data){
          chart['abscissaValue'] = JSON.parse(data[0]).map(Number);
          // Input the downloaded results into the chart object
          const downloadedOrdinates = data.slice(1);
          // Convert the results into javascript objects/variables
          for(let i=0; i < downloadedOrdinates.length; i++) {
            // Convert each String into a var/object usable by JS
            downloadedOrdinates[i] = JSON.parse(downloadedOrdinates[i]);
          }
          // Add the downloaded data with the right format into the chart object
          chart['selectedResultsValues'] = downloadedOrdinates;

          /* Resolve the promise returning the chart object with the downloaded
          data into it */
          resolve(chart);
        }).catch(err => {reject(err)});
      });
    },
    // Return true if the object contains one-dimension variables
    containsOneDimensionVariables: function(array) {
      return this.containsNDimensionsVariables(array, 1);
    },
    // Return true if the object contains two-dimensions variables
    containsTwoDimensionVariables: function(array) {
      return this.containsNDimensionsVariables(array, 2);
    },
    containsNDimensionsVariables: function(array, n) {
      for(let i=0; i<array.length; i++) {
        if(this.getNumberOfDimensions(array[i]) == n) return true;
      }
      return false;
    },
    // Return the number of mathematical dimensions of a given variable
    getNumberOfDimensions: function(variable) {
      let nbOfDimensions = 0;
      let temp = variable;

      if(temp.length === undefined) {
        nbOfDimensions = 1;
      } else {
        while (typeof(temp) === 'object') {
          nbOfDimensions++;
          temp = temp[0];
        }
      }

      return nbOfDimensions;
    },
    /* Generate an info message for the user in case she chose to display
    one dimension variables and matrices on the same chart */
    generateInfoMsgForMultipleNbDimensionsChoice: function(chart) {
      const oneDimensionVariables = [];
      const matrices = [];
      /* Sort the variables, the one dimension on one side, multiple
      dimensions on the other side */
      for (let i=0;i<chart.selectedResultsValues.length;i++) {
        const nbDimensions = this.getNumberOfDimensions(chart.selectedResultsValues[i]);
        if(nbDimensions > 1) matrices.push(chart.selectedResults[i]);
        else oneDimensionVariables.push(chart.selectedResults[i]);
      }

      // Construct an info text to let the user know which data to choose
      let info = "Error, impossible to generate a graph with one dimension variables and matrices.\nOne dimension variables are : ";
      for (var index in oneDimensionVariables) {
        if (oneDimensionVariables.hasOwnProperty(index)) {
          info = `${info}${oneDimensionVariables[index]}, `
        }
      }
      info = info.slice(0, info.length -2);
      info = info + '\nMatrices are : ';

      for (var index in matrices) {
        if (matrices.hasOwnProperty(index)) {
          info = `${info}${matrices[index]}, `
        }
      }
      info = info.slice(0, info.length -2);
      info = info +'.';

      return info;
    }
  }
}
</script>

<style lang="scss">
.generatedChart {
  margin-bottom: 20px;
}
</style>
