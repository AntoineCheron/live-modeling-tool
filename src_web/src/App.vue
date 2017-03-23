<template>
  <div id="app">
    <nav class="navbar navbar-default">
      <div class="container-fluid">
        <div class="navbar-header">
          <a class="navbar-brand" href="#">Live modeling tool for OSUR - Boussinesq Simulation</a>
        </div>
      </div>
    </nav>
    <input-files-section></input-files-section>
    <charts-section></charts-section>
    <results-selection></results-selection>
  </div>
</template>

<script>
import InputFilesSection from './Inputs.vue'
import ChartsSection from './Charts.vue'
import ResultsSelection from './Results.vue'

export default {
  name: 'app',
  data () {
    return {
    }
  },
  components : {
    InputFilesSection,
    ChartsSection,
    ResultsSelection
  },
}
</script>

<style lang="scss">

</style>
