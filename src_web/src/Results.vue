<template>
  <section id="resultsSelection">
    <div class="col-lg-3">
      <div class="col-lg-12">
        <h1>Results list</h1>
        <ul id="resultsList">
          <!-- Results downloaded from the API will be placed here -->
        </ul>
      </div>
      <div class="col-lg-12">
        <h1>Charts selection</h1>
        <!-- TODO : create this dynamic component that let the user add and remove charts on the go -->
      </div>
    </div>
  </section>
</template>

<script>
export default {
  name: 'results',
  data () {
    return {
      msg: 'this is the results view',
    }
  }
}
</script>

<style lang="scss"></style>
