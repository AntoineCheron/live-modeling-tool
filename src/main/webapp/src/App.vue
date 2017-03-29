<template>
  <div id="app">
    <nav class="navbar navbar-default">
      <div class="container-fluid">
        <div class="navbar-header">
          <a class="navbar-brand" href="#">Live modeling tool for OSUR - Boussinesq Simulation</a>
        </div>
      </div>
    </nav>
    <input-files-section class="col-lg-3"></input-files-section>
    <results class="col-lg-9"></results>
  </div>
</template>

<script>
import InputFilesSection from './Inputs.vue'
import Results from './Results.vue'

export default {
  name: 'app',
  data () {
    return {
    }
  },
  components : {
    InputFilesSection,
    Results
  },
}
</script>

<style lang="scss">

</style>
