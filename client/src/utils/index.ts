import { extent } from "d3-array";
import { scaleLinear } from "d3-scale";

export function normalizeArray(
  list: number | number[],
  range: [number, number],
  domain?: [number | undefined, number | undefined]
): number | number[] {
  if (!Array.isArray(list)) {
    return list;
  }
  if (!domain) {
    domain = extent(list);
  }
  if (typeof domain[0] !== "number" || typeof domain[1] !== "number") {
    return list;
  }
  if (domain[0] === domain[1] && domain[0] >= range[0] && domain[0] <= range[1]) {
    return list;
  }
  var scale = scaleLinear()
    .domain(domain as [number, number])
    .range(range);
  return list.map((n) => scale(n));
}

export function toNumber(n: number | string | null | undefined, def: number = 0): number {
  if (n == null) {
    return def;
  }
  const v = parseFloat(n.toString());
  return isNaN(v) ? def : v;
}

export function toExtent<T>(
  data: T[],
  getter: (d: T) => [number, number] | [undefined, undefined] | undefined
): [number, number] | [undefined, undefined] {
  return data.reduce<[number, number] | [undefined, undefined]>(
    (acc, d) => {
      const current = getter(d);
      if (current) {
        if (acc[0] == null || (current[0] != null && current[0]! < acc[0]!)) {
          acc[0] = current[0];
        }
        if (acc[1] == null || (current[1] != null && current[1]! > acc[1]!)) {
          acc[1] = current[1];
        }
      }
      return acc;
    },
    [undefined, undefined]
  );
}

export async function downloadFile(data: any, name: string) {
  const json = JSON.stringify(data);
  const blob = new Blob([json], { type: "application/json" });
  const href = await URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = href;
  link.download = name + ".json";
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}
