export const DEFAULT_COLORWAY = [
  "#4e79a7",
  "#f28e2c",
  "#e15759",
  "#76b7b2",
  "#59a14f",
  "#edc949",
  "#af7aa1",
  "#ff9da7",
  "#9c755f",
  "#bab0ab",
];

export const ARRAY_DISTANCE_METRICS: { [key: string]: (value: number, values: number[]) => number } = {
  euclidean_difference: (value, values) => Math.sqrt(values.reduce((acc, cur) => acc + Math.pow(cur - value, 2), 0)),
  mean: (value, values) => values.reduce((acc, cur) => acc + Math.abs(cur), 0) / values.length,
  mean_rel_diff: (value, values) => values.reduce((acc, cur) => acc + Math.abs(cur), 0) / values.length / value,
  mean_rel_diff_log10: (value, values) =>
    values.reduce((acc, cur) => acc + Math.log10(Math.abs(cur)), 0) / values.length / Math.log10(value),
  max: (value, values) => Math.max(...values),
  max_difference: (value, values) => Math.abs(Math.abs(value) - Math.abs(Math.max(...values))),
  min: (value, values) => Math.min(...values),
  min_difference: (value, values) => Math.abs(Math.abs(value) - Math.abs(Math.min(...values))),
};

export const isProxySymbol = Symbol("isProxy");

export const DEFAULT_PLOTLY_COLORSCALE = {
  colorscale: [
    ["0", "#D0EBFF"],
    ["0.5", "#49a7bd"],
    ["1.0", "#31347a"],
  ],
  reversescale: false,
};
