<template>
  <section id="results">
    <h1>Results</h1>
    <button @click="runSimulation" class="btn btn-default">Run simulation</button>
    <div class="loader" v-if="simulationRunning"></div>
    <charts :charts="chartsToDisplay"></charts>
    <charts-selection
    :charts="charts"
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
  data() {
    return {
      charts: [],
      chartsToDisplay: [],
      nextChartId: 1,
      simulationRunning: false
    }
  },
  methods: {
    runSimulation: function() {
      WS.getSimulate()
      .then(this.simulationEnded)
      .catch(function(err){console.error(err.statusText)});
      this.simulationRunning = true;
    },
    simulationEnded: function() {
      this.simulationRunning = false;
    },
    addResultToChart: function(chart, result) {
      const index = _.indexOf(this.charts, chart);
      this.charts[index].selectedResults.push(result);
    },
    removeResultFromChart: function(chart, result) {
      const index = _.indexOf(this.charts, chart);
      this.charts[index].selectedResults = _.without(this.charts[index].selected_results, result);
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
      this.charts.push({id:this.nextChartId, type: 'line chart', selectedResults: [], abscissa: '', displayed: false});
      this.nextChartId = this.nextChartId+1;
    },
    generateChart: function(chart) {
      // Change the displayed value of chart in this.charts to true
      const i = this.charts.indexOf(chart);
      console.log(i);
      chart.displayed = true;
      this.charts[i] = chart;
      console.log(this.charts);
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
    }
  },
  components: {
    Charts,
    ChartsSelection
  }
}
</script>

<style lang="scss">
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
