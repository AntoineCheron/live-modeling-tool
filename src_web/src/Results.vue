<template>
  <section id="results">
    <h1>Results</h1>
    <button v-on:click="runSimulation">Run simulation</button>
    <div class="loader" v-if="simulationRunning"></div>
    <charts></charts>
    <charts-selection></charts-selection>
  </section>
</template>

<script>
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
      WS.getSimulate(this.simulationEnded);
      this.simulationRunning = true;
    },
    simulationEnded: function() {
      this.simulationRunning = false;
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
