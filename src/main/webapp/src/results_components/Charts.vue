<template>
  <section id="charts">
    <div class="col-lg-8">
      <h1>Charts</h1>
      <div v-for="chart in charts">
        <individual-chart :chart="chart" @error="removeChart(chart)"></individual-chart>
      </div>
    </div>
  </section>
</template>

<script>
import IndividualChart from "./IndividualChart.vue"
import _ from "underscore"

export default {
  name: 'charts',
  props: {
    'charts': {
      type: Array,
      required: true
    }
  },
  data () {
    return {
    }
  },
  methods: {
    removeChart: function(chart) {
      this.$emit('errorDisplayingChart', chart);
    },
  },
  components: {
    IndividualChart
  }
}
</script>

<style lang="scss">
</style>
