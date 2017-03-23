<template>
  <section id="resultsSelection">
    <div class="col-lg-3">
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
        <chart-selection v-for="chart in charts" v-bind:chart-data="resultsFormat"
        v-bind:chart-list="chartTypes" v-bind:chartType="chart.type"
        v-on:selectedResult="result => {addResultToChart(chart, result)}"
        v-on:unselectedResult="result => {removeResultFromChart(chart, result)}" v-on:remove="removeChart(chart)"
        v-on:selectedChartType="type => {setChartType(chart, type)}" v-on:generate="generateChart(chart)"></chart-selection>
      </div>
      <button class="btn btn-default" id="addChartButton" v-on:click="addChart">Add a new chart</button>
    </div>
  </section>
</template>

<script>
import _ from 'underscore'
import ChartSelection from './results_components/ChartComponent.vue'

export default {
  name: 'results',
  data () {
    return {
      resultsFormat: [
        'Q',
        'S',
        'QS',
      ],
      charts: [
        {id:1, type: 'default', selectedResults: []},
        {id:2, type: 'default', selectedResults: []},
      ],
      chartTypes: ['pie chart', 'flow chart', 'line chart'],
      nextChartId: 3,
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
    },
    setChartType: function(chart, type) {
      // Set the var type of the chart
      const index = _.indexOf(this.charts, chart);
      this.charts[index].type = type;
    },
    generateChart: function(chart) {
      this.$emit('generateChart', chart);
    },
    addChart: function() {
      this.charts.push({id:this.nextChartId, type: 'default', selectedResults: []});
      this.nextChartId = this.nextChartId+1;
      console.log(this.charts);
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
</style>
