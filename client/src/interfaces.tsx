export interface IParticle {
  selected?: boolean;
  collection: string;
  index?: number;
  structure: string;
  images?: string[];
  embedding?: number[];
  properties?: {
    clusters?: string | null | undefined;
    [key: string]: string | number | boolean | null | undefined;
  };
  projection: { [key: string]: number[] };
  plotData?: {
    color?: string;
  };
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

export interface ICollection<T extends IParticle = IParticle> {
  data: T[];
  name: string;
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
