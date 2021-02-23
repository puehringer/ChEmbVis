import * as React from "react";
import { Figure } from "react-plotly.js";
import { ICollection, IParticle, IParticleSelection, IPlotOptions } from "../interfaces";
import { PlotComponent } from "./PlotComponent";
import lodashGet from "lodash.get";
import { extent } from "d3-array";
import { normalizeArray, toExtent, toNumber } from "../utils";
import merge from "lodash.merge";

export interface IProjectionPlotProps {
  title: string;
  collections: ICollection[];
  xAccessor: string;
  yAccessor: string;
  options: IPlotOptions;
  setHover(particle: IParticle | null): void;
  selected: IParticleSelection;
  setSelected(selected: IParticleSelection): void;
  children?: React.ReactNode;
}

export const ScatterPlot = React.memo(
  ({
    title,
    collections,
    xAccessor,
    yAccessor,
    options,
    setHover,
    selected,
    setSelected,
    children,
  }: IProjectionPlotProps) => {
    const [figureState, setFigureState] = React.useState<Figure | null>(null);
    const [innerHover, setInnerHover] = React.useState<IParticle | null>(null);

    React.useEffect(() => {
      const timeout = setTimeout(() => {
        setFigureState((figureState) => {
          const annotatedCollections = collections.map(({ data, name, plotOptions }) => {
            const colorBy = plotOptions?.colorBy || options.colorBy;
            const opacityBy = plotOptions?.opacityBy || options.opacityBy;
            const sizeBy = plotOptions?.sizeBy || options.sizeBy;

            const color = colorBy ? data.map((p) => toNumber(p.properties?.[colorBy])) : undefined;
            const constantOpacity = typeof opacityBy === "number" ? opacityBy : undefined;
            const opacity = opacityBy
              ? data.map((p) => constantOpacity || toNumber(p.properties?.[opacityBy]))
              : undefined;
            // opacity: opacityBy
            // ? normalizeArray(
            //     points.map((p) => p.properties?.[opacityBy] as number),
            //     0.3,
            //     0.7
            //   )
            // : defaultOpacity,
            const size = sizeBy ? data.map((p) => toNumber(p.properties?.[sizeBy], 5)) : undefined;
            // TODO: Move outside
            // @ts-ignore
            const colorExtent: [number, number] = color ? extent(color) : undefined;
            // @ts-ignore
            const opacityExtent: [number, number] = opacity ? extent(opacity) : undefined;
            // @ts-ignore
            const sizeExtent: [number, number] = size ? extent(size) : undefined;

            return {
              data,
              name,
              plotOptions,
              colorExtent,
              color,
              opacity,
              opacityExtent,
              size,
              sizeExtent,
            };
          });

          const colorExtent = toExtent(annotatedCollections, (d) => d.colorExtent);
          const opacityExtent = toExtent(annotatedCollections, (d) => d.opacityExtent);
          const sizeref = options.sizeBy
            ? (2.0 * Math.max(...annotatedCollections.map(({ sizeExtent }) => sizeExtent?.[1] || 0))) / 5 ** 2
            : undefined;

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
            data: annotatedCollections.map(({ data, name, color, opacity, size, plotOptions }, i) => {
              const existingData = (figureState?.data.find((d) => d.name === name) as Partial<Plotly.PlotData>) || {};

              const groupBy = plotOptions?.groupBy || options.groupBy;
              const colorBy = plotOptions?.colorBy || options.colorBy;
              const groups = groupBy ? data.map((p) => p.properties?.[groupBy]?.toString() || "N/A") : undefined;

              return {
                ...existingData,
                x: data.map((p) => lodashGet(p, xAccessor)),
                y: data.map((p) => lodashGet(p, yAccessor)),
                text: data.map((p) => p.structure),
                hoverinfo: "none",
                name,
                type: "scattergl",
                mode: "markers",
                marker: {
                  ...(existingData.marker || {}),
                  color,
                  cmin: colorExtent?.[0],
                  cmax: colorExtent?.[1],
                  // opacity: opacity ? normalizeArray(opacity, [0.0, 0.7], opacityExtent) : 0.5,
                  opacity: opacity || 0.5,
                  symbol: i,
                  size: size ?? 5,
                  sizeref,
                  sizemin: 2,
                  sizemax: 5,
                  colorbar: colorBy
                    ? {
                        title: colorBy.length > 10 ? colorBy.slice(0, 7) + "..." : colorBy,
                        titleside: "top",
                        thickness: 10,
                        outlinewidth: 0,
                        lenmode: "fraction",
                        len: 0.75,
                        yanchor: "middle",
                      }
                    : undefined,
                  colorscale: "Portland",
                },
                transforms: (groupBy
                  ? [
                      {
                        type: "groupby",
                        groups,
                        styles: Array.from(new Set(groups)).map((target, i) => ({
                          target,
                          value: { marker: { symbol: i } },
                        })),
                      },
                    ]
                  : []) as Plotly.Transform[],
              };
            }),
          };
          return newFigureState;
        });
      }, 300);

      return () => {
        clearTimeout(timeout);
      };
    }, [title, collections, xAccessor, yAccessor, options]);

    React.useEffect(() => {
      setFigureState((figureState) =>
        figureState
          ? {
              ...figureState,
              data: figureState.data.map((data) => ({
                ...data,
                selectedpoints: selected?.[data.name!]?.map((s) => s.index!),
              })),
            }
          : null
      );
    }, [selected]);

    React.useEffect(() => {
      const timeout = setTimeout(() => setHover(innerHover), 50);

      return () => {
        clearTimeout(timeout);
      };
    }, [innerHover, setHover]);

    const getPointsFromEvent = (e: Readonly<Plotly.PlotSelectionEvent> | null): IParticleSelection => {
      const points = e?.points || [];
      return points.length > 0
        ? points.reduce<{ [key: string]: IParticle[] }>((acc, p) => {
            // @ts-ignore
            const color = (p?.fullData as Plotly.Data)?.marker?.color;
            const collection = collections?.[p.curveNumber];
            const particle: IParticle | null = collection?.data?.[p.pointIndex];
            if (particle) {
              particle.plotData = {
                ...particle.plotData,
                color: typeof color === "string" ? color : undefined,
              };
              if (!acc[collection.name]) {
                acc[collection.name] = [];
              }
              acc[collection.name].push(particle);
            }
            return acc;
          }, {})
        : null;
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
              setSelected(null);
            }}
            onClick={(e) => {
              setSelected(getPointsFromEvent(e));
              setInnerHover(null);
            }}
            onHover={(e) => {
              setInnerHover(Object.values(getPointsFromEvent(e)!).flat()?.[0]);
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
