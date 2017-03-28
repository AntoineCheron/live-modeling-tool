<template>
  <section id="chartsSelection">
    <div class="col-lg-4">
      <div class="col-lg-12">
        <h1>Results list</h1>
        <ul id="resultsList">
          <li v-for="result in resultsFormat">
            {{ result }}
          </li>
        </ul>
      </div>
      <div class="col-lg-12">
        <h1>Charts selection</h1>
        <chart-selection v-for="chart in charts"
        :chart-data="resultsFormat"
        :chart-list="chartTypes"
        :chart-object="chart"
        @selectedResult="result => {addResultToChart(chart, result)}"
        @unselectedResult="result => {removeResultFromChart(chart, result)}"
        @selectedAbscissa="result => {setAbscissa(chart, result)}"
        @remove="removeChart(chart)"
        @selectedChartType="type => {setChartType(chart, type)}"
        @generate="generateChart(chart)"
        @update="updateChart(chart)">
        </chart-selection>
      </div>
      <button class="btn btn-default" id="addChartButton" @click="addChart">Add a new chart</button>
    </div>
  </section>
</template>

<script>
import WS from '../WebServices/WebServices.js'
import _ from 'underscore'
import ChartSelection from './ChartSelectionComponent.vue'

export default {
  name: 'results',
  created () {
    WS.getResulsFormat()
    .then(this.updateResultsFormat)
    .catch(function(err){console.error(err.statusText)});
  },
  data () {
    return {
      resultsFormat: [],
      charts: [],
      chartTypes: ['pie chart', 'flow chart', 'line chart'],
      nextChartId: 1,
    }
  },
  methods: {
    addResultToChart: function(chart, result) {
      const index = _.indexOf(this.charts, chart);
      this.charts[index].selectedResults.push(result);
    },
    removeResultFromChart: function(chart, result) {
      const index = _.indexOf(this.charts, chart);
      this.charts[index].selectedResults = _.without(this.charts[index].selected_results, result);
    },
    removeChart: function(chart) {
      this.charts = _.without(this.charts, chart);
      this.$emit('removeChart', chart);
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
    generateChart: function(chart) {
      this.$emit('generateChart', chart);
    },
    updateChart: function(chart) {
      this.$emit('updateChart', chart);
    },
    addChart: function() {
      this.charts.push({id:this.nextChartId, type: 'line chart', selectedResults: [], abscissa: '', displayed: false});
      this.nextChartId = this.nextChartId+1;
    },
    // Callback for the getRequest method
    updateResultsFormat: function(responseText) {
      // Verify that the result is not empty before modifying the value
      const newResultsFormat = JSON.parse(responseText);
      if (newResultsFormat !== undefined) {
        this.resultsFormat = newResultsFormat;
      }
    }
  },
  components: {
    ChartSelection
  }
}
</script>

<style lang="scss">
  #addChartButton{
    margin-top: 20px;
  }

  form {
    display: inline-block;
    margin-right: 8px;
    margin-left: 4px;
  }

  select {
    display: block;
    margin-top: 10px;
    margin-bottom: 10px;
  }

  .chartElement {
    margin-bottom: 20px;
  }
</style>
