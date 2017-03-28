<template>
  <section id="results">
    <h1>Results</h1>
    <button @click="runSimulation">Run simulation</button>
    <div class="loader" v-if="simulationRunning"></div>
    <charts :charts="charts"></charts>
    <charts-selection
    @removeChart="chart => {removeChart(chart)}"
    @generateChart="chart => {addChart(chart)}"
    @updateChart="chart => {updateChart(chart)}"></charts-selection>
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
    removeChart: function(chart) {
      // Retrieve the index of chart into this.charts
      let index = null;
      for(let i=0;i<this.charts.length;i++) {
        if(this.charts[i].id === chart.id) {
          index = i;
          break;
        }
      }

      // Remove the chart at the index
      if(index !== null){
        this.charts = _.without(this.charts, this.charts[index]);
      }
    },
    addChart: function(chart) {
      chart.displayed = true;
      this.charts.push(chart);
    },
    updateChart: function(chart) {
      // Retrieve the index of chart into this.charts
      let index = null;
      for(c in this.charts) {
        if(c.id === chart.id) index = this.charts.indexOf(c);
      }

      if (index !== null){
        // Replace the chart at index in the this.charts Array
        this.charts[index] = chart

        // It will automatically update the view
      }
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
