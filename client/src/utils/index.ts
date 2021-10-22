import { extent } from "d3-array";
import { scaleLinear } from "d3-scale";

export function normalizeArray<T extends number | number[]>(
  list: T,
  range: [number, number],
  domain?: [number | undefined, number | undefined]
): T {
  if (typeof list === "number" || !Array.isArray(list)) {
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
  return list.map((n) => scale(n)) as any;
}

export function toNumber(n: number | string | null | boolean | undefined, def: number = 0): number {
  if (n == null) {
    return def;
  }
  if (n === true) {
    return 1;
  } else if (n === false) {
    return 0;
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

export function stringifyWithoutCycles(data: any): string {
  return JSON.stringify(data, (key, value) => {
    if (key === "knn_particles") {
      return null;
    }
    return value;
  });
}

export async function downloadFile(blob: Blob, name: string) {
  const href = await URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = href;
  link.download = name;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

export async function downloadJSONFile(data: any, name: string) {
  const json = stringifyWithoutCycles(data);
  const blob = new Blob([json], { type: "application/json" });
  downloadFile(blob, name + ".json");
}

export async function downloadCSVFile(data: string, name: string) {
  const blob = new Blob([data], { type: "text/csv" });
  downloadFile(blob, name + ".csv");
}

/**
 * Injects a cached getter function into an object.
 * @param object Object to be inejected into.
 * @param key Key of the newly injected getter.
 * @param valueCreator Function returning the value of the injected getter.
 * @param options Options.
 */
export function injectGetter<T>(object: object, key: PropertyKey, valueCreator: () => T, {
  cache = true
}: {
  cache?: boolean;
} = {}): void {
  Object.defineProperty(object, key, {
    get: function () {
      const value = valueCreator();
      if(cache) {
        // Remove the getter function.
        delete (object as any)[key];
        // Assign the actual value and return it.
        // return ((object as any)[key] = value);
        // Inject the getter again instead, now without the "replace" logic.
        injectGetter(object, key, () => value, { cache: false });
      }
      return value;
    },
    configurable: true,
    enumerable: true,
  });
}
