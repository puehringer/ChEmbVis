import * as React from "react";
import { Figure } from "react-plotly.js";
import { IParticle } from "../interfaces";
import { PlotComponent } from "./PlotComponent";
import lodashGet from "lodash.get";
import memoizeOne from "memoize-one";

function _getSelectedIndices<T>(data: T[], selection: T[]): number[] {
  const lookup = data.reduce<Map<T, number>>((acc, cur, i) => {
    acc.set(cur, i);
    return acc;
  }, new Map<T, number>());
  return selection.map((s) => lookup.get(s)!).filter((s) => s != null);
}
const getSelectedIndices = memoizeOne(_getSelectedIndices);

export interface IProjectionPlotProps {
  title?: string;
  marker?: Partial<Plotly.PlotMarker>;
  markerFunction?: (points: IParticle[], index: number) => Partial<Plotly.PlotMarker>;
  transformsFunction?: (points: IParticle[], index: number) => Plotly.DataTransform[];
  additionalTraces?: Plotly.Data | Plotly.Data[] | null;
  additionalTracesFunction?: (
    points: IParticle[],
    get: (p: IParticle | null, axis: "x" | "y") => number
  ) => Plotly.Data[];
  transforms?: Plotly.DataTransform[];
  xAccessor: string;
  yAccessor: string;
  particles: IParticle[];
  additionalParticles: IParticle[];
  setHover(particle: IParticle | null): void;
  selected: IParticle[];
  setSelected(selected: IParticle[]): void;
  children?: React.ReactNode;
}

export const ProjectionPlot = React.memo(
  ({
    xAccessor,
    yAccessor,
    title,
    particles,
    additionalParticles,
    markerFunction,
    transformsFunction,
    setHover,
    selected,
    setSelected,
    additionalTraces,
    additionalTracesFunction,
    children,
  }: IProjectionPlotProps) => {
    const [figureState, setFigureState] = React.useState<Figure | null>(null);
    const [innerHover, setInnerHover] = React.useState<IParticle | null>(null);

    const transforms = React.useMemo(() => transformsFunction?.(particles, 0), [particles, transformsFunction]);
    const marker = React.useMemo(() => markerFunction?.(particles, 0), [particles, markerFunction]);
    const additionalTransforms = React.useMemo(() => transformsFunction?.(particles, 1), [
      particles,
      transformsFunction,
    ]);
    const additionalMarker = React.useMemo(() => markerFunction?.(particles, 1), [particles, markerFunction]);
    const additionalComputedTraces = React.useMemo(
      () => additionalTracesFunction?.(particles, (p, axis) => lodashGet(p, axis === "x" ? xAccessor : yAccessor)),
      [particles, xAccessor, yAccessor, additionalTracesFunction]
    );

    const mainTrace = React.useMemo<Plotly.Data>(() => {
      return {
        // ...(figureState.data?.[0] || {}) as Partial<PlotData>,
        x: particles.map((p) => lodashGet(p, xAccessor)),
        y: particles.map((p) => lodashGet(p, yAccessor)),
        text: particles.map((p) => p.structure),
        hoverinfo: "none",
        name: "Embedding",
        type: "scattergl",
        mode: "markers",
        marker,
        transforms,
      };
    }, [particles, xAccessor, yAccessor, marker, transforms]);

    const mainTraceWithSelection = React.useMemo<Plotly.Data>(() => {
      const selectedpoints = getSelectedIndices(particles, selected);

      return {
        ...mainTrace,
        selectedpoints: selectedpoints.length > 0 ? selectedpoints : undefined,
      };
    }, [mainTrace, selected, particles]);

    const secondTrace = React.useMemo<Plotly.Data>(() => {
      return {
        x: additionalParticles.map((p) => lodashGet(p, xAccessor)),
        y: additionalParticles.map((p) => lodashGet(p, yAccessor)),
        text: additionalParticles.map((p) => p.structure),
        hoverinfo: "none",
        name: "Interpolated",
        type: "scatter",
        mode: "markers",
        marker: additionalMarker,
        transforms: additionalTransforms,
      };
    }, [additionalParticles, xAccessor, yAccessor, additionalMarker, additionalTransforms]);

    const secondTraceWithSelection = React.useMemo<Plotly.Data>(() => {
      const additionalSelectedpoints = getSelectedIndices(additionalParticles, selected);

      return {
        ...secondTrace,
        selectedpoints: additionalSelectedpoints.length > 0 ? additionalSelectedpoints : undefined,
      };
    }, [secondTrace, selected, additionalParticles]);

    React.useEffect(() => {
      setFigureState((figureState) => {
        const newFigureState: Figure = {
          frames: [],
          layout: {
            ...(figureState?.layout || {}),
            dragmode: "lasso",
            hovermode: "closest",
            autosize: true,
            legend: {
              x: 1,
              xanchor: "right",
              y: 1,
            },
            title,
            margin: {
              // l: 0,
              r: 0,
              // b: 0,
              // t: 25,
              pad: 4,
            },
          },
          data: [
            mainTraceWithSelection,
            secondTraceWithSelection,
            ...(additionalTraces ? (Array.isArray(additionalTraces) ? additionalTraces : [additionalTraces]) : []),
            ...(additionalComputedTraces || []),
          ],
        };

        // return merge(figureState, newFigureState);

        return newFigureState;
      });
    }, [mainTraceWithSelection, secondTraceWithSelection, additionalTraces, additionalComputedTraces, title]);

    React.useEffect(() => {
      const timeout = setTimeout(() => setHover(innerHover), 50);

      return () => {
        clearTimeout(timeout);
      };
    }, [innerHover, setHover]);

    const getPointsFromEvent = (e: Readonly<Plotly.PlotSelectionEvent> | null): IParticle[] => {
      return (
        e?.points
          .map((p) => {
            // @ts-ignore
            const color = (p?.fullData as Plotly.Data)?.marker?.color;
            let particle: IParticle | null = null;
            if (p.curveNumber === 0) {
              particle = particles?.[p.pointIndex];
            } else if (p.curveNumber === 1) {
              particle = additionalParticles?.[p.pointIndex];
            } else if (p.curveNumber === 2) {
              const structure = p.data.text?.[p.pointIndex];
              if (structure) {
                // TODO: Make sure a particle with structure is ok...
                // @ts-ignore
                particle = {
                  structure,
                };
              }
            }
            if (particle) {
              particle.plotData = {
                ...particle.plotData,
                color: typeof color === "string" ? color : undefined,
              };
            }
            // I don't know why typescript does not check the .filter(Boolean)
            return particle!;
          })
          .filter(Boolean) || []
      );
    };

    return (
      <>
        {figureState ? (
          <PlotComponent
            style={{
              width: "100%",
              height: "100%",
            }}
            data={figureState.data}
            layout={figureState.layout}
            config={{
              displaylogo: false,
              responsive: true,
              showLink: false,
              // showEditInChartStudio: true,
              // plotlyServerURL: "https://chart-studio.plotly.com"
            }}
            onSelected={(e) => {
              setSelected(getPointsFromEvent(e));
            }}
            onDeselect={() => {
              setSelected([]);
            }}
            onClick={(e) => {
              setSelected(getPointsFromEvent(e));
              setInnerHover(null);
            }}
            onHover={(e) => {
              setInnerHover(getPointsFromEvent(e)[0]);
            }}
            onUnhover={() => {
              setInnerHover(null);
            }}
            onUpdate={(figure) => {
              // Inline save the figure to save the zoom
              figureState.data = figure.data;
              figureState.frames = figure.frames;
              figureState.layout = figure.layout;
            }}
          />
        ) : null}
        {children}
      </>
    );
  }
);
