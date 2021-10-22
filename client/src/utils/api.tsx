import { stringifyWithoutCycles } from ".";
import { IInterpolatedParticle, IParticle, IRegistry, IServerCollection } from "../interfaces";

function fetchRaw({ url, data, method = "POST" }: { url: string; data?: any; method?: string }): Promise<Response> {
  return fetch(url, {
    headers: {
      "Content-Type": "application/json",
    },
    // @ts-ignore
    // mode: '*cors', // no-cors, *cors, same-origin
    method,
    redirect: "follow",
    body: stringifyWithoutCycles(data),
  }).then(async (res) => {
    if (!res.ok) {
      throw Error((await res.json().catch(() => null))?.message || res.statusText);
    }
    return res;
  });
}

function fetchJSON<T>(url: string, data: any): Promise<T> {
  return fetchRaw({ url, data }).then((res) => res.json());
}

function fetchText(url: string, data: any): Promise<string> {
  return fetchRaw({ url, data }).then((res) => res.text());
}

export function runMSO(options: {
  structure: string;
  iterations?: number;
  num_swarms?: number;
  num_part?: number;
  v_min?: number;
  v_max?: number;
  inertia_weight?: number;
  phi1?: number;
  phi2?: number;
  phi3?: number;
  objectives: any[];
}): Promise<IServerCollection<IParticle>> {
  return fetchJSON("/api/pso/", options);
}

export function embedStructures(data: {
  structures: string[];
  include_embedding?: boolean;
}): Promise<IServerCollection<IParticle>> {
  return fetchJSON<IServerCollection<IParticle>>("/api/embedding/", data);
}

export function getChemblUMAPEmbedding(): Promise<IParticle[]> {
  return fetchRaw({
    url: "/api/embedding/",
    method: "GET",
  }).then((res) => res.json());
}

export function interpolateStructures(
  structures: string[],
  maxSamples?: number
): Promise<IServerCollection<IInterpolatedParticle>> {
  return fetchJSON("/api/interpolation/", {
    structures,
    maxSamples,
  });
}

export function getImageURL(
  structure: string,
  substructure: string | null = null,
  align: string | null = null
): string {
  return `/api/image/?structure=${encodeURIComponent(structure)}${
    substructure ? `&substructure=${encodeURIComponent(substructure)}` : ""
  }${align ? `&align=${encodeURIComponent(align)}` : ""}`;
}

export function getReducedImages(
  structures: string[],
  method: "single" | "murcko" | "mcs" | "similarity" | "auto" = "auto"
): Promise<string | null> {
  return fetchText("/api/image/", {
    structures,
    method,
  }).catch(() => null);
}

export function computeProjectionsWithModels(
  particles: Pick<IParticle, 'structure' | 'embedding'>[],
  models: { [key: string]: string }
): Promise<{
  projections: IServerCollection["projections"];
  projection: { [key: string]: number[][] };
  additional?: { [key: string]: number[][] | null };
}> {
  return fetchJSON("/api/projection/models", {
    particles: particles.map((p) => ({structure: p.structure, embedding: p.embedding})),
    models,
  });
}

export function hasSubstructureMatch(
  structures: string[],
  smarts: string
): Promise<{
  validity: { [key: string]: boolean };
  counts: { [key: string]: number };
}> {
  return fetchJSON("/api/mol/substructures/", {
    structures,
    smarts,
  });
}

export function getTanimotoSimilarity(
  structures: string[],
  reference: string,
  fingerprint: string
): Promise<{
  tanimoto: { [key: string]: number };
}> {
  return fetchJSON("/api/mol/tanimoto/", {
    structures,
    reference,
    fingerprint,
  });
}

export function getNeighborSamples(
  structure: string,
  samples: number,
  method: string,
  scale: number
): Promise<IServerCollection<IInterpolatedParticle>> {
  return fetchJSON("/api/sampling/", {
    structure,
    samples,
    method,
    scale,
  });
}

export function getMatchedMolecularPairs(options: {
  structure: string;
  min_variable_size?: number;
  max_variable_size?: number;
  min_constant_size?: number;
  min_radius?: number;
  min_pairs?: number;
  substructure?: string;
}): Promise<{
  structures: string[];
}> {
  return fetchJSON("/api/mmp/", options);
}

export function getStonedSelfies(options: {
  structure: string;
  substructure: string;
  random_samples?: number;
  max_mutations?: number;
}): Promise<{
  structures: string[];
}> {
  return fetchJSON("/api/stoned_selfies/", options);
}

export function getRegistry(): Promise<IRegistry> {
  return fetchRaw({
    url: "/api/registry/",
    method: "GET",
  }).then((res) => res.json());
}
