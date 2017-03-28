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
        @generate="generateChart(chart)">
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
  props: {
    'charts': {
      type: Array,
      required: true
    },
  },
  data () {
    return {
      resultsFormat: [],
      charts: [],
      chartTypes: ['pie chart', 'flow chart', 'line chart'],
    }
  },
  methods: {
    addResultToChart: function(chart, result) {
      this.$emit('addResultToChart', chart, result);
    },
    removeResultFromChart: function(chart, result) {
      this.$emit('removeResultFromChart', chart, result);
    },
    setChartType: function(chart, type) {
      this.$emit('setChartType', chart, type);
    },
    setAbscissa: function(chart, abscissa) {
      this.$emit('setAbscissa', chart, abscissa);
    },
    generateChart: function(chart) {
      this.$emit('generateChart', chart);
    },
    addChart: function() {
      this.$emit('addChart');
    },
    removeChart: function(chart) {
      this.$emit('removeChart', chart);
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
