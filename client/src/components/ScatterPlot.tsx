import * as React from "react";
import { Figure } from "react-plotly.js";
import { ICollection, IParticle, IParticleSelection, IPlotOptions } from "../interfaces";
import { PlotComponent } from "./PlotComponent";
import lodashGet from "lodash.get";
import { extent } from "d3-array";
import { normalizeArray, toExtent, toNumber } from "../utils";
import groupBy from "lodash.groupby";
import { scaleLinear, scaleSymlog } from "d3-scale";
import { color } from "d3-color";

const TRAJECTORY_TRACE_NAME = "Trajectories";

export interface IProjectionPlotProps {
  title: string;
  collections: ICollection[];
  xAccessor: string;
  yAccessor: string;
  options: IPlotOptions;
  hover: IParticle | null;
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
    hover: debouncedHover,
    setHover,
    selected,
    setSelected,
    children,
  }: IProjectionPlotProps) => {
    const [figureState, setFigureState] = React.useState<Figure | null>(null);
    const [innerHover, setInnerHover] = React.useState<IParticle | null>(null);

    // Use the most up-to-date hover if possible.
    const hover = innerHover || debouncedHover;

    React.useEffect(() => {
      const traces: Partial<Plotly.PlotData>[] = ((() => {
        // void clusters; // Line only exists to avoid eslint-disable-next-line react-hooks/exhaustive-deps
        // return null;
        const connectBy = options.connectBy;
        if (!connectBy || connectBy.length === 0 || (hover == null && !selected)) {
          return null;
        }

        const filteredCollections: Map<ICollection | undefined, IParticle[]> = new Map();
        if (selected) {
          Object.entries(selected).forEach(([name, selection]) => {
            const collection = collections.find((c) => c.name === name);
            filteredCollections.set(collection, selection);
          });
        }

        filteredCollections.delete(undefined);
        if (!hover && filteredCollections.size === 0) {
          return null;
        }

        return [
          (hover ? [collections.find((c) => c.data.includes(hover)), [hover]] : [undefined, undefined]) as [
            ICollection | undefined,
            IParticle[] | undefined
          ],
          ...Array.from(filteredCollections.entries()),
        ].map(([c, selection], i) => {
          if (!c || !selection) {
            return null;
          }
          const particles = c!.data;
          const isHover = i === 0;

          const selectedInstancesByConnectBy = connectBy.map(
            (c) => new Set(selection.map((p) => p?.properties?.[c]?.toString()).filter((id) => id != null))
          );

          const filteredParticles = particles.filter((p) =>
            connectBy.every((connectBy, i) =>
              selectedInstancesByConnectBy[i].has(p.properties?.[connectBy]?.toString())
            )
          );

          if (filteredParticles.length === 0) {
            return null;
          }

          // TOOD: This could be memoized.
          const groups = Object.entries(
            groupBy(filteredParticles, (p) =>
              connectBy.map((connectBy) => `${connectBy}:${p.properties?.[connectBy]}`).join(", ")
            )
          );

          const allInstances: (IParticle | null)[] = [];
          const allSizes: (number | null)[] = [];
          const allColors: (string | null)[] = [];

          const lineOpacityScaling = scaleSymlog().range([0.1, 0.8]).domain([particles.length, 0]);
          const hoverColor = color(/* hover?.plotData?.color || */ "rgb(0,0,0)")?.darker()!;

          const lineColor = hoverColor.copy();
          lineColor.opacity = lineOpacityScaling(filteredParticles.length);
          const markerBorderColor = "gray";
          hoverColor.opacity = 0.5;
          const sizeScaling = scaleLinear().range([6, 12]);
          // Cool plotly optimization: instead of creating many traces for lines, create a single trace with NaN separators.
          // See https://www.somesolvedproblems.com/2019/05/how-to-make-plotly-faster-with-many.html
          for (let [key, instances] of groups) {
            instances = instances.filter((p) => lodashGet(p, xAccessor) != null && lodashGet(p, yAccessor) != null);
            const instanceScaler = sizeScaling.domain([0, instances.length]);
            allInstances.push(...instances);
            allInstances.push(null);
            allSizes.push(...instances.map((_, i) => instanceScaler(i)));
            allSizes.push(null);
            allColors.push(...instances.map((p, i) => (p === hover || p.selected ? "gold" : "lightgray")));
            allColors.push("darkblue");
          }

          return {
            type: isHover ? "scatter" : "scattergl",
            mode: "lines+markers",
            x: allInstances.map((p) => lodashGet(p, xAccessor) ?? NaN),
            y: allInstances.map((p) => lodashGet(p, yAccessor) ?? NaN),
            name: TRAJECTORY_TRACE_NAME,
            hoverinfo: "skip",
            opacity: hover ? (isHover ? 1.0 : 0.05) : 0.2,
            marker: {
              color: allColors,
              size: allSizes,
              line: {
                color: markerBorderColor.toString(),
                width: 2,
              },
              opacity: 1,
            },
            line: {
              color: lineColor.toString(),
              width: 2,
              shape: "spline",
              // opacity: 0.5
            },
            showlegend: false,
          } as Partial<Plotly.PlotData>;
        });
      })()?.filter(Boolean) || []) as Partial<Plotly.PlotData>[];

      setFigureState((figureState) => {
        if (!figureState) {
          return figureState;
        }
        const existing = figureState.data.find((d) => d.name === TRAJECTORY_TRACE_NAME);
        if (traces.length === 0 && !existing) {
          return figureState;
        }
        const data = [
          ...figureState.data
            .filter((d) => d.name !== TRAJECTORY_TRACE_NAME)
            .map((d) => ({ ...d, opacity: traces.length > 0 ? 0.2 : undefined })),
          ...traces,
        ];
        return { ...figureState, data };
      });
      // const result: Plotly.Data[] = [];
      // groups.forEach(([group, instances]) => {
      //   // const c = color(
      //   //   instances.find((i) => i.plotData?.color)?.plotData?.color || "#000000"
      //   // )!;
      //   const c = color(hover?.plotData?.color || "#000000")!;
      //   const sizeScaling = scaleLinear()
      //     .domain([0, instances.length])
      //     .range([8, 16]);
      //   c.opacity = hover && instances.includes(hover) ? 0.5 : 0.3;
      //   if (true || instances.length > 20 || groups.length > 20) {
      //     result.push({
      //       type: "scattergl",
      //       mode: "lines+markers",
      //       x: instances.map((p) => get(p, "x")),
      //       y: instances.map((p) => get(p, "y")),
      //       name: group,
      //       hoverinfo: "all",
      //       marker: {
      //         color: "lightblue",
      //         // @ts-ignore
      //         size: instances.map((p, i) => sizeScaling(i)),
      //         line: {
      //           color: c.toString(),
      //           width: 2,
      //         },
      //         opacity: 1,
      //         // line: 'green'
      //       },
      //       line: {
      //         color: c.toString(),
      //         // shape: "spline",
      //       },
      //       showlegend: false,
      //     });
      //   } else {
      //     // const widthScaling = scaleLinear()
      //     //   .domain([0, instances.length])
      //     //   .range([1, 6]);
      //     // instances.slice(1).forEach((instance, i) => {
      //     //   const instancePair = [instances[i], instance];
      //     //   result.push({
      //     //     type: "scattergl",
      //     //     mode: "lines+markers",
      //     //     x: instancePair.map((p) => get(p, "x")),
      //     //     y: instancePair.map((p) => get(p, "y")),
      //     //     name: group,
      //     //     hoverinfo: "all",
      //     //     marker: {
      //     //       color: 'red',
      //     //       size: 10
      //     //     },
      //     //     line: {
      //     //       width: widthScaling(i),
      //     //       color: c.toString(),
      //     //       shape: "spline",
      //     //       // dash: 'dot',
      //     //     },
      //     //     // legendgroup: `${group}`,
      //     //     // showlegend: i === 0,
      //     //     showlegend: false,
      //     //   });
      //     // });
      //   }
      // });

      // return result;
    }, [collections, xAccessor, yAccessor, options.connectBy, hover, selected]);

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
                  opacity: opacity ? normalizeArray(opacity, [0.1, 0.9]/* , opacityExtent */) : 0.5,
                  // opacity: opacity || 0.5,
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
                        styles: Array.from(new Set(groups))
                          .sort((a, b) => a.localeCompare(b))
                          .map((target, i) => ({
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
