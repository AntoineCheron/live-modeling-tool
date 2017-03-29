<template>
  <section id="results">
    <h1>Results</h1>
    <div class="horizontalContainer">
      <button @click="runSimulation" v-if="!simulationRunning" class="btn btn-default">Run simulation</button>
      <button @click="runSimulation" v-if="simulationRunning" class="btn btn-default" disabled>Simulation running</button>
      <div class="loader" v-if="simulationRunning"></div>
      <h3 v-if="simulationRunning">Simulation running for <span class="blue">{{simulationTimer.minutes}}</span>min<span class="orange">{{simulationTimer.seconds}}</span>seconds</h3>
    </div>
    <charts :charts="chartsToDisplay"></charts>
    <charts-selection
    :charts="charts"
    :resultsFormat="resultsFormat"
    @addResultToChart="(chart, result) => {addResultToChart(chart, result)}"
    @removeResultFromChart="(chart, result) => {removeResultFromChart(chart, result)}"
    @setChartType="(chart, type) => {setChartType(chart, type)}"
    @setAbscissa="(chart, abscissa) => {setAbscissa(chart, abscissa)}"
    @addChart="addChart"
    @removeChart="chart => {removeChart(chart)}"
    @generateChart="chart => {generateChart(chart)}"></charts-selection>
  </section>
</template>

<script>
import _ from "underscore"
import WS from "./WebServices/WebServices.js"
import Charts from "./results_components/Charts.vue"
import ChartsSelection from "./results_components/ChartsSelection.vue"

export default {
  created () {
    this.updateResultsFormat();
  },
  data() {
    return {
      charts: [],
      chartsToDisplay: [],
      nextChartId: 1,
      simulationRunning: false,
      simulationTimer: {
        minutes: 0,
        seconds: 0
      },
      simulationTimeInterval: null,
      resultsFormat: [],
    }
  },
  methods: {
    runSimulation: function() {
      // Request the server to launch the simulation
      WS.getSimulate()
      .then(this.simulationEnded)
      .catch(function(err){
        this.simulationEnded();
        console.error(err.statusText)
      });
      // Turns the simulation running var to true to modify the view
      this.simulationRunning = true;
      // Launch the timer
      const Parent = this;
      Parent.simulationTimer.seconds = 0;
      Parent.simulationTimer.minutes = 0;
      this.simulationTimeInterval = setInterval(function(){
        Parent.simulationTimer.seconds++;
        if(Parent.simulationTimer.seconds == 60) {
          Parent.simulationTimer.minutes++;
          Parent.simulationTimer.seconds = 0;
        }
      }, 1000);
    },
    simulationEnded: function() {
      this.simulationRunning = false;
      clearInterval(this.simulationTimeInterval);
      this.updateResultsFormat();
    },
    addResultToChart: function(chart, result) {
      const index = _.indexOf(this.charts, chart);
      this.charts[index].selectedResults.push(result);
    },
    removeResultFromChart: function(chart, result) {
      const index = _.indexOf(this.charts, chart);
      this.charts[index].selectedResults = _.without(this.charts[index].selectedResults, result);
    },
    setChartType: function(chart, type) {
      // Set the var type of the chart
      const index = _.indexOf(this.charts, chart);
      this.charts[index].type = type;
    },
    setAbscissa: function(chart, abscissa) {
      // Set the var type of the chart
      const index = _.indexOf(this.charts, chart);
      this.charts[index].abscissa = abscissa;
    },
    addChart: function() {
      this.charts.push({id:this.nextChartId, type: 'default', selectedResults: [], abscissa: '', displayed: false});
      this.nextChartId = this.nextChartId+1;
    },
    generateChart: function(chart) {
      // Change the displayed value of chart in this.charts to true
      const i = this.charts.indexOf(chart);
      chart.displayed = true;
      this.charts[i] = chart;
      // Add the chart object into the chartsToDisplay array
      this.chartsToDisplay.push(chart);
    },
    removeChart: function(chart) {
      // First : remove the chart in the chartsToDisplay array
      // Retrieve the index of chart into this.chartsToDisplay
      let index = null;
      for(let i=0;i<this.chartsToDisplay.length;i++) {
        if(this.chartsToDisplay[i].id === chart.id) {
          index = i;
          break;
        }
      }
      // Remove the chart at the index in chartsToDisplay
      if(index !== null){
        this.chartsToDisplay = _.without(this.chartsToDisplay, this.chartsToDisplay[index]);
      }

      /* Second : remove the chart in the charts array. It is easier because
      chart is synchronised with charts, but not with chartsToDisplay*/
      this.charts = _.without(this.charts, chart);
    },
    // Callback for the getRequest method
    updateResultsFormat: function() {
      const ParentContext = this;
      WS.getResultsFormat()
      .then(function(responseText) {
        // Verify that the result is not empty before modifying the value
        const newResultsFormat = JSON.parse(responseText);
        if (newResultsFormat !== undefined) {
          ParentContext.resultsFormat = newResultsFormat;
        }
      })
      .catch(function(err){console.error(err.statusText)});
    }
  },
  components: {
    Charts,
    ChartsSelection
  }
}
</script>

<style lang="scss">
.horizontalContainer {
  display: flex;
  align-items: center;
}

.blue {
  color: #268bd2;
}

.orange {
  color: #cb4b16;
}

.loader {
  display: inline-block;
  border: 6px solid #f3f3f3; /* Light grey */
  border-top: 6px solid #3498db;
  border-radius: 50%;
  width: 30px;
  height: 30px;
  animation: spin 2s linear infinite;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% {transform: rotate(360deg); }
}
</style>
