<template>
  <div class="chartElement">
    <p>Please select the results you want to see on the chart : </p>

    <form v-for="result in chartData">
      <input type="checkbox" v-bind:id="result" v-bind:value="result"
      v-bind:checked="chartObject.selectedResults.indexOf(result) != -1"
      v-on:change="changeOnResult($event.target, result)">
      <label v-bind:for="result">{{ result }}</label>
    </form>
    <br />
    <select v-bind:value="chartObject.type" v-on:input="selectedChartType($event.target.options[$event.target.selectedIndex].value)">
      <option value="default">Please select the type of the chart</option>
      <option v-for="type in chartList" v-bind:value="type">{{ type }}</option>
    </select>

    <button class="btn btn-success" v-on:click="generate">Generate chart</button>
    <button class="btn btn-default" v-on:click="remove">Remove</button>
  </div>
</template>

<script>
export default {
  props : {
    'chartData': {
      type: Array,
      required: true
    },
    'chartList': {
      type: Array,
      required: true
    },
    'chartObject': {
      type: Object,
      required: true
    },
  },
  methods: {
    changeOnResult: function(eventTarget, result) {
      eventTarget.checked ? this.$emit('selectedResult', result) : this.$emit('unselectedResult', result);
    }, remove: function() {
      this.$emit('remove');
    }, selectedChartType: function(type){
      this.$emit('selectedChartType', type);
    }, generate: function() {
      this.$emit('generate');
    }
  }
}
</script>

<style lang="scss" scoped>

</style>
