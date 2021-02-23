export interface IParticle {
  selected?: boolean;
  collection: string;
  index?: number;
  structure: string;
  embedding?: number[];
  properties?: {
    clusters?: string | null | undefined;
    [key: string]: string | number | null | undefined;
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

export interface IPlotOptions {
  colorBy: string | null;
  opacityBy: number | string | null;
  groupBy: string | null;
  connectByValues: string[] | null;
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

export declare type IParticleSelection = { [key: string]: IParticle[] } | null;

export enum EActiveTabs {
  EMBEDDING,
  INTERPOLATION,
}
