<template>
  <div class="chartElement">
    <p>Please select the results you want to see on the chart : </p>

    <form v-for="result in chartData">
      <input type="checkbox" :value="result"
      :checked="chartObject.selectedResults.indexOf(result) != -1"
      @change="changeOnResult($event.target, result)">
      <label :for="result">{{ result }}</label>
    </form>
    </br>
    <select :value="chartObject.abscissa" @input="selectedAbscissa($event.target.options[$event.target.selectedIndex].value)">
      <option value="">Please select the abscissa</option>
      <option v-for="result in chartData" :value="result">{{ result }}</option>
    </select>
    <select :value="chartObject.type" @input="selectedChartType($event.target.options[$event.target.selectedIndex].value)">
      <option value="default">Please select the type of the chart</option>
      <option v-for="type in chartList" :value="type">{{ type }}</option>
    </select>

    <button class="btn btn-success" @click="generate" v-if="!chartObject.displayed">Generate chart</button>
    <button class="btn btn-default" @click="remove">Remove</button>
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
    }, selectedAbscissa: function(abscissa){
      this.$emit('selectedAbscissa', abscissa);
    }, generate: function() {
      this.$emit('generate');
    }
  }
}
</script>

<style lang="scss" scoped>
  select {
    margin-top: 10px;
  }
</style>
