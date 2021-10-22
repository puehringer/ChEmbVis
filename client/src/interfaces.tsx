export interface INearestNeighbors {
  distance_metric: string;
  knn_dist: number[];
  knn_ind: number[];
  knn_particles: IParticle[];
}

export interface IParticle {
  selected?: boolean;
  collection: string;
  index?: number;
  structure: string;
  original_structure?: string;
  images?: string[];
  embedding?: { [key: string]: number[] };
  nearest_neighbors?: {
    [key: string | symbol]: INearestNeighbors;
  };
  clusters?: {
    [key: string | symbol]: {
      distance_metric: string;
      label: number;
    };
  };
  properties: {
    [key: string | symbol]: string | number | boolean | null | undefined;
  };
  projection: { [key: string]: number[] };
}

export function isParticle(p: any): p is IParticle {
  return p.structure && p.projection;
}

export const DEFAULT_COLLECTION = "Particles";
export const DEFAULT_CHEMBL_COLLECTION = "ChEMBL";

export interface IPlotOptions {
  colorBy: string | null;
  opacityBy: number | string | null;
  groupBy: string | null;
  connectBy: string[] | null;
  sizeBy: string | null;
}

export interface IEnabledProjection {
  label: string;
  value: string;
  projection: string | null;
  plotOptions: Partial<IPlotOptions>;
}

export interface IServerCollection<T extends IParticle = IParticle> {
  data: T[];
  projections?: {
    [key: string]: {
      trustworthiness?: number;
      trustworthiness_additional?: number;
      explained_variance?: number;
      model?: string;
    };
  };
}

export interface ICollection<T extends IParticle = IParticle> extends IServerCollection<T> {
  name: string;
  type?: "neighborhoodSampling";
  hidden?: boolean;
  plotOptions?: Partial<IPlotOptions>;
}

export interface IInterpolatedParticle extends IParticle {
  scaffold: boolean;
}

export interface IObjective {
  description: string;
  name: string;
  weight: number;
  desirability: { x: number; y: number }[];
  additional_args?: { [key: string]: string };
}
export interface IRegistry {
  objectives: IObjective[];
}

export declare type IParticleSelection = { [key: string]: IParticle[] } | null;

export enum EActiveTabs {
  EMBEDDING,
  INTERPOLATION,
}
