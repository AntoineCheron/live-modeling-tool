<template>
  <section id="charts">
    <div class="col-lg-8">
      <h1>Charts</h1>
      <div class="generatedChart" v-for="chart in chartsWithOptions">
        <highcharts :options="chart.options"></highcharts>
      </div>
    </div>
  </section>
</template>

<script>
import WS from "../WebServices/WebServices.js"
import Vue from "vue"
import VueHighcharts from "vue-highcharts"
import AsyncComputed from "vue-async-computed"
import _ from "underscore"

Vue.use(VueHighcharts);
Vue.use(AsyncComputed);

export default {
  name: 'charts',
  props: {
    'charts': {
      type: Array,
      required: true
    }
  },
  asyncComputed: {
    chartsWithOptions: function() {
      const parent = this;
      return new Promise(function(resolve, reject) {
        // Verify that charts isn't empty. If so return an empty Array
        if(parent.charts.length === 0) resolve([]);

        const charts = parent.charts;
        const computedCharts = [];
        /* Because this.generateOptionObjectFromChart is a promise and there
        potentially are several charts, we need to use a PromiseArray
        to store all the promise in order to execute a Promise.All then.*/
        const promiseArray = [];

        // Call this.generateOptionObjectFromChart for each chart to display
        for (let i=0; i < charts.length ; i++) {
          promiseArray.push(parent.generateOptionObjectFromChart(charts[i]));
        }

        /* After all options objects have been computed, put them all into the
        computedCharts array that will be used to render the view */
        // Retrieve the charts object to use it into the promise.all

        Promise.all(promiseArray).then(function(optionsArray){
          /* Add the options object into each chart and then the chart into the
          computedCharts array */
          for(let i=0; i < optionsArray.length; i++) {
            const newChart = _.clone(charts[i]);
            newChart['options'] = optionsArray[i];
            computedCharts.push(newChart);
          }

          // Finish :D, return the computedCharts array
          resolve(computedCharts);
        }).catch(err => {reject(err)});
      });
    }
  },
  data () {
    return {
    }
  },
  methods: {
    generateOptionObjectFromChart: function(chart) {
      /* This promise is used to return the options object that will be
      built after having downloaded all the needed elements */
      return new Promise(function (resolve, reject) {
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
        have been downloaded. Thus, we can create the options object used by
        Highcharts */
        Promise.all(Downloads).then(function(data){
          const xAxis = JSON.parse(data[0]).map(Number);

          // Prepare the yAxis and title that will be input into the options object
          // The results (ordinates) are all the data array minus the first element
          const downloadedOrdinates = data.slice(1);
          const yAxisSeries = [];
          let maxSerieSize = 0;
          let yAxisTitle = "";

          for (let i = 0; i < downloadedOrdinates.length ; i++) {
            const serieLength = JSON.parse(downloadedOrdinates[i]).length;
            if( serieLength > maxSerieSize){
              maxSerieSize = serieLength;
            }
            // create an yAxis serie correctly formatted for each data
            const tempYAxisSerie = {
              name: chart.selectedResults[i],
              data: JSON.parse(downloadedOrdinates[i]).map(Number)
            };
            // Add it into the yAxisSeries array
            yAxisSeries.push(tempYAxisSerie);
            // Continue to compute the title
            yAxisTitle += `${chart.selectedResults[i]}, `;
          }
          // Make title looking nicer, removing the , and space at the end
          yAxisTitle = yAxisTitle.slice(0, yAxisTitle.length -2);

          /* Important step : HERE WE HAVE ALL THE ELEMENTS TO BUILD
          THE OPTIONS OBJECT, LET'S DO THIS*/
          const options = {
            title: {
              text: `Chart ${chart.id}`
            },
            xAxis: {
              categories: xAxis
            },
            yAxis: {
              title: {
                text: yAxisTitle
              },
              plotLines: [{
                value: 0,
                width: 1,
                color: '#808080'
              }]
            },
            legend: {
              layout: 'vertical',
              align : 'right',
              verticalAlign: 'middle',
              borderWidth: 0
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
          resolve(options);

        }).catch(function(err){reject(err)});

      });
    }
  }
}
</script>

<style lang="scss">
.generatedChart {
  margin-bottom: 20px;
}
</style>
