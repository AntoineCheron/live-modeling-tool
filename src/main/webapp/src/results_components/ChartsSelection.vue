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
  props: {
    'charts': {
      type: Array,
      required: true
    },
    'resultsFormat': {
      type: Array
    }
  },
  data () {
    return {
      charts: [],
      chartTypes: ['pie', 'line', 'area', 'areaspline', 'bar', 'column', 'scatter', 'spline'],
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
