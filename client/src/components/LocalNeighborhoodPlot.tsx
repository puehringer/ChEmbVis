import { buildRanking, NumberColumn, Ranking, Taggle } from "lineupjs";
import * as React from "react";
import Select from "react-select/creatable";
import { ICollection, INearestNeighbors, IParticle, IParticleSelection } from "../interfaces";
import { normalizeArray, toNumber } from "../utils";
import { PlotComponent, PLOTLY_CONFIG } from "./PlotComponent";
import { LineupWrapper } from "./ranking/LineupWrapper";
import { StructureImageColumn } from "./ranking/StructureImageColumn";
import { StructureImage } from "./StructureImage";

function getKNNByMetric(
  particles: IParticle[],
  reference: IParticle,
  getter: (particle: IParticle) => number
): INearestNeighbors {
  const values = particles.map((particle) => getter(particle));
  const referenceValue = getter(reference);
  const difference = values
    .map((value, index) => ({ diff: Math.abs(value - referenceValue), index }))
    .sort((a, b) => a.diff - b.diff)
    .slice(0, 50);
  return {
    distance_metric: "absolute_difference",
    knn_ind: difference.map(({ index }) => index),
    knn_dist: difference.map(({ diff }) => diff),
    knn_particles: difference.map(({ index }) => particles[index]),
  };
}

export function LocalNeighborhoodPlot({
  selected,
  collection,
}: //   setHover,
//   hover,
{
  selected: IParticle | null | undefined;
  collection: ICollection | null | undefined;
  setHover(hover: IParticle | null): void;
  hover: IParticle | null;
}) {
  const [hover, setHover] = React.useState<IParticle | null>(null);
  const [selection, setSelection] = React.useState<Set<IParticle>>(new Set());
  const [referenceEmbedding, setReferenceEmbedding] = React.useState<string>("");
  const [enabledEmbeddings, setEnabledEmbeddings] = React.useState<string[]>([]);

  React.useEffect(() => {
    setHover(null);
    setSelection(new Set());
  }, [selected]);

  const allAvailableProperties = React.useMemo(() => Object.keys(selected?.properties || {}), [selected]);
  /* const allAvailableEmbeddings = React.useMemo(
    () =>
      Object.entries({
        ...(selected?.nearest_neighbors || {}),
        // TODO: Enable to allow properties as "distance"
        // ...Object.fromEntries(
        //   collection && selected
        //     ? allAvailableProperties
        //         .filter(([key, value]) => typeof value === "number")
        //         .map(([key, value]) => [
        //           key,
        //           getKNNByMetric(collection.data, selected, (p) => p.properties[key] as number),
        //         ])
        //     : []
        // ),
      }),
    [selected, collection, allAvailableProperties]
  ); */

  const filteredAvailableEmbeddings = React.useMemo(
    () =>
      enabledEmbeddings.length === 0 || !selected
        ? []
        : enabledEmbeddings
            .map((e) => [e, selected.nearest_neighbors?.[e]!] as [string, INearestNeighbors])
            .filter(([key, value]) => key && value),
    [enabledEmbeddings, selected]
  );

  const getRankingBuilders = React.useCallback(() => {
    // For each available embedding, create a ranking
    return filteredAvailableEmbeddings.map(([key, value], i) => {
      const rankingBuilder = buildRanking();
      // rankingBuilder.column('_rank');
      rankingBuilder.column("structure");
      if (i === 0) {
        rankingBuilder.column("Occurances");
      }
      rankingBuilder.column(`Distance_${key}`);
      rankingBuilder.sortBy(`Distance_${key}`);
      // rankingBuilder.allColumns();
      return rankingBuilder;
    });
  }, [filteredAvailableEmbeddings]);

  const adjustRankings = React.useCallback((lineup: Taggle, rankings: Ranking[]) => {
    rankings.forEach((ranking) => {
      const distanceColumn = ranking.flatColumns.find(
        (col) => col instanceof NumberColumn && col.desc.label.toLowerCase().startsWith("distance")
      ) as NumberColumn | undefined;
      // distanceColumn
      distanceColumn?.setFilter({
        filterMissing: true,
        min: -Infinity,
        max: Infinity,
      });
      // distanceColumn?.on("filterChanged", (prev, current) => console.log(current));
    });
  }, []);

  // Compute the "histogram" of index occurances --> the more often it occurs, the "better" the neighbor is
  const indexOccurances = React.useMemo(
    () =>
      filteredAvailableEmbeddings.reduce<{ [index: number]: number }>((acc, [key, value]) => {
        value.knn_ind.forEach((v) => {
          acc[v] = (acc[v] || 0) + 1;
        });
        return acc;
      }, {}),
    [filteredAvailableEmbeddings]
  );

  const selectionCollections = React.useMemo(
    () => [
      {
        name: "Selection",
        data: (selection.size > 0
          ? Array.from(selection)
          : Array.from(new Set(filteredAvailableEmbeddings.map(([key, value]) => value.knn_particles).flat()))
        ).map((p) => ({
          ...p,
          // Additionally store the original particle as reference for selection synchronization
          originalParticle: p,
          properties: {
            Occurances: filteredAvailableEmbeddings
              .filter(([key, value]) => value.knn_particles.includes(p))
              .map(([key, value]) => key!) as any,
            ...p.properties,
            ...filteredAvailableEmbeddings.reduce<{ [key: string]: number | null }>((acc, [key, value]) => {
              const index = value.knn_particles.indexOf(p);
              acc[`Distance_${key}`] = index >= 0 ? value.knn_dist[index] : null;
              return acc;
            }, {}),
          },
        })),
      },
    ],
    [selection, filteredAvailableEmbeddings]
  );

  const setSelectedFromLineup = React.useCallback(
    (selection: IParticleSelection) =>
      setHover(
        // @ts-ignore
        Object.values(selection || {}).flat()?.[0]?.originalParticle
      ),
    []
  );

  const normalizedDistances = React.useMemo(
    () =>
      Object.fromEntries(
        filteredAvailableEmbeddings.map(([key, value]) => [key, normalizeArray(value.knn_dist, [0, 1])])
      ),
    [filteredAvailableEmbeddings]
  );
  const referenceIndexToDist = React.useMemo(() => {
    const referenceNN = selected?.nearest_neighbors?.[referenceEmbedding];
    return referenceNN && normalizedDistances[referenceEmbedding]
      ? Object.fromEntries(referenceNN.knn_ind.map((cur, i) => [cur, normalizedDistances[referenceEmbedding][i]]))
      : null;
  }, [selected, referenceEmbedding, normalizedDistances]);

  const plotDataLines = React.useMemo(() => {
    const data: Partial<Plotly.PlotData>[] = [];
    // For all embedding pairs (y-values)
    filteredAvailableEmbeddings.forEach(([key, value], i) => {
      if (i < 1) {
        return;
      }
      const [previousKey, previousValue] = filteredAvailableEmbeddings[i - 1];
      const previousInd = new Set(previousValue.knn_ind);
      // console.log(selection.size > 0
      //   ? value.knn_particles.map((v, i) => (selection.has(v) ? i : null)).filter((i) => i != null)
      //   : undefined);
      // Go through all nearest neighbors
      value.knn_ind.forEach((ind, i) => {
        // Check if the nearest neighbor is also in the previous track
        if (previousInd.has(ind) && (selection.size === 0 || selection.has(value.knn_particles[i]))) {
          // Retrieve the x values for the current and previous track
          const x = normalizedDistances[key][i];
          const previousX = normalizedDistances[previousKey][previousValue.knn_ind.indexOf(ind)];
          // Add a line to the plot for exactly this pair
          data.push({
            x: [previousX, x],
            y: [previousKey, key],
            type: "scatter",
            mode: "lines",
            showlegend: false,
            legendgroup: "Pairwise Connections",
            opacity: Math.abs(previousX - x) / (value.knn_ind.length / 5),
            line: {
              color: "gray",
            },
          });
        }
      });
    });
    return data;
  }, [filteredAvailableEmbeddings, normalizedDistances, selection]);

  const plotDataSelection = React.useMemo<Partial<Plotly.PlotData>[]>(() => {
    const selectedXValues = hover
      ? filteredAvailableEmbeddings.map(([key, value], i) => {
          const index = value.knn_particles.indexOf(hover);
          return index >= 0 ? normalizedDistances[key][index] : NaN;
        })
      : null;

    if (selectedXValues) {
      return [
        {
          x: selectedXValues.map((v) => (isNaN(v) ? 1 : v)),
          y: filteredAvailableEmbeddings.map(([key, value]) => key),
          type: "scatter",
          mode: "lines+markers",
          showlegend: false,
          line: {
            color: "gold",
          },
          marker: {
            color: selectedXValues.map((v) => (isNaN(v) ? "red" : "gold")),
            // size: selectedXValues.map((v) => isNaN(v) ? 0 : 10)
          },
        },
      ];
    }
    return [];
  }, [filteredAvailableEmbeddings, hover, normalizedDistances]);

  const plotDataScatter = React.useMemo(() => {
    return filteredAvailableEmbeddings.map(([key, value], i) => {
      // By default, color each point by the number of occurances in the different embeddings
      let color: Plotly.Color = value.knn_ind.map((ind) => indexOccurances[ind] / filteredAvailableEmbeddings.length);
      //   If we have a reference, color each point by the position in the reference embedding
      if (referenceIndexToDist) {
        color = value.knn_ind.map((ind) => referenceIndexToDist[ind] || NaN);
      } else if (referenceEmbedding) {
        color = value.knn_particles.map((p) => toNumber(p.properties[referenceEmbedding]));
      }

      return {
        x: normalizedDistances[key],
        y: value.knn_dist.map(() => key),
        text: value.knn_dist.map((dist) => `${dist} (${value.distance_metric})`),
        hoverinfo: "none" as const,
        name: key,
        type: "scatter" as const,
        mode: "markers" as const,
        customdata: value.knn_particles.map((p) => p.index!),
        selectedpoints:
          selection.size > 0
            ? value.knn_particles.map((v, i) => (selection.has(v) ? i : null)).filter((i) => i != null)
            : undefined,
        showlegend: false,
        showscale: true,
        // line: {
        //   width: 0
        // },
        marker: {
          //   opacity: value.knn_ind.map((ind) => indexOccurances[ind] / filteredAvailableEmbeddings.length),
          opacity: 0.4,
          //   color: value.knn_ind.map((ind) => indexOccurances[ind] / filteredAvailableEmbeddings.length),
          color,
          coloraxis: "coloraxis",
        },
      };
    });
  }, [
    filteredAvailableEmbeddings,
    indexOccurances,
    normalizedDistances,
    referenceEmbedding,
    referenceIndexToDist,
    selection,
  ]);

  const data: Partial<Plotly.PlotData>[] = [...plotDataLines, ...plotDataSelection, ...plotDataScatter];

  const getPointsFromEvent = (e: Pick<Plotly.PlotMouseEvent, "points">): IParticle[] => {
    console.log(e?.points);
    return collection
      ? Array.from(
          new Set(e?.points?.sort((a, b) => +a.x! - +b.x!)?.map((p) => collection.data[p.customdata as number]))
        ).filter(Boolean)
      : [];
  };

  return (
    <div
      className="d-flex flex-column"
      style={{
        position: "relative",
        // overflowY: 'scroll',
        // minHeight: 500
      }}
    >
      <p className="text-muted">
        Inspect the nearest neighbors in the embedding spaces of the selected or hovered structure. Each embedding shows
        the top-N nearest neighbors according to its distance metric. Generally, the color encodes in how many
        embeddings a specific neighbor is also listed as neighbor, or if a reference embedding is given, it shows the
        color of the position in the reference embedding (or black if it is no common neighbor).
      </p>
      {selected && collection ? (
        <>
          <div className="row">
            <div className="col d-flex flex-column align-items-center justify-content-center">
              <strong>Selected</strong>
              {selected?.structure ? <StructureImage style={{ width: "100%" }} structure={selected.structure} align={undefined /* selected?.structure */} /> : null}
            </div>
            <div className="col d-flex flex-column align-items-center">
              <strong>Hovered</strong>
              {hover?.structure ? (
                <StructureImage style={{ width: "100%" }} structure={hover.structure} align={undefined /* selected?.structure */} />
              ) : null}
            </div>
            <div className="col d-flex flex-column align-items-center">
              <strong>Similarity Map</strong>
              {selected?.structure && hover?.structure ? (
                <StructureImage
                  style={{ width: "100%" }}
                  structure={[selected.structure, hover.structure]}
                  align={selected?.structure}
                />
              ) : null}
            </div>
          </div>
          <div className="row">
            <div className="col-md-8 mb-3">
              <label>Enabled embeddings</label>
              <Select<{ label: string; value: string }, true>
                menuPosition="fixed"
                isMulti
                value={enabledEmbeddings.map((p) => ({
                  label: p,
                  value: p,
                }))}
                onChange={(e) => {
                  setEnabledEmbeddings(e.map((p) => p.value));
                }}
                options={[
                  ...Object.entries(selected?.nearest_neighbors || {}).map(([key, value]) => ({
                    label: `${key} (${value.distance_metric})`,
                    value: key,
                  })),
                  ...allAvailableProperties.map((property) => ({
                    label: `Property: ${property}`,
                    value: property,
                  })),
                ]}
                closeMenuOnSelect={false}
              />
            </div>
            <div className="col-md-4 mb-3">
              <label>Color by</label>
              <Select<{ label: string; value: string }, false>
                menuPosition="fixed"
                value={
                  referenceEmbedding
                    ? {
                        label: referenceEmbedding,
                        value: referenceEmbedding,
                      }
                    : null
                }
                isClearable={true}
                onChange={(e) => {
                  setReferenceEmbedding(e?.value || "");
                }}
                options={[
                  ...Object.entries(selected?.nearest_neighbors || {})
                    .map(([key, value]) => key)
                    .map((key) => ({
                      label: key,
                      value: key,
                    })),
                  ...allAvailableProperties.map((property) => ({
                    label: `Property: ${property}`,
                    value: `property=${property}`,
                  })),
                ]}
              />
            </div>
          </div>
          {data.length > 0 ? (
            <>
              <PlotComponent
                className="mt-3 mb-5"
                style={{ width: "100%", height: filteredAvailableEmbeddings.length * 40 + 50 }}
                data={data}
                layout={{
                  dragmode: "lasso",
                  hovermode: "closest",
                  autosize: true,
                  legend: {
                    // x: 1,
                    // xanchor: "right",
                    // y: 1,
                    orientation: "h",
                  },
                  // @ts-ignore
                  coloraxis: {
                    // cmin: referenceIndexToDist ? undefined : 0,
                    // cmax: referenceIndexToDist ? undefined : 1,
                    colorscale: "Portland",
                    // reversescale: true,
                    colorbar: {
                      title: "Distance",
                      titleside: "top",
                      thickness: 10,
                      outlinewidth: 0,
                      lenmode: "fraction",
                      len: 0.75,
                      yanchor: "center",
                    },
                  },
                  scene: {
                    aspectmode: "data",
                  },
                  margin: {
                    // l: 250,
                    r: 0,
                    b: 0,
                    t: 50,
                    // pad: 4,
                  },
                  yaxis: {
                    categoryorder: "array",
                    categoryarray: filteredAvailableEmbeddings.map(([key, value]) => key).reverse(),
                  },
                }}
                config={PLOTLY_CONFIG}
                onSelected={(e) => {
                  setSelection(new Set(getPointsFromEvent(e)));
                }}
                onHover={(e) => {
                  setHover(getPointsFromEvent(e)?.[0]);
                }}
                onUnhover={() => {
                  setHover(null);
                }}
              />
              <div style={{ flex: "1 1 500px", display: "flex", overflow: "auto" }}>
                {selectionCollections ? (
                  <LineupWrapper
                    collections={selectionCollections}
                    setSelected={setSelectedFromLineup}
                    // align={hover?.structure}
                    getRankingBuilders={enabledEmbeddings.length <= 1 ? undefined : getRankingBuilders}
                    adjustRankings={enabledEmbeddings.length === 0 ? undefined : adjustRankings}
                  />
                ) : null}
              </div>
            </>
          ) : null}
        </>
      ) : null}
    </div>
  );
}
