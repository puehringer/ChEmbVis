import * as React from "react";
import { Modal, Button } from "react-bootstrap";
import { StructureCard } from "./components/StructureCard";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { StructureImage } from "./components/StructureImage";
import { ICollection, IParticle, IParticleSelection } from "./interfaces";
import groupByLodash from "lodash.groupby";
import { embedStructures, getTanimotoSimilarity } from "./utils/api";
import { PlotComponent, PLOTLY_CONFIG } from "./components/PlotComponent";

export function getClustersFromParticle(particle: IParticle): string[] {
  return particle?.properties?.clusters?.split(";").filter(Boolean) || [];
}

export const ClusterSidePanel = React.memo(
  ({
    clusters,
    setClusters,
    selected,
    setSelected,
  }: {
    clusters: { [key: string]: IParticle[] };
    setClusters(name: string, clusters: IParticle[] | null): void;
    selected: IParticleSelection;
    setSelected(selected: IParticleSelection): void;
  }) => {
    const selectedParticles = React.useMemo<IParticle[]>(
      () => (selected ? ([] as IParticle[]).concat(...Object.values(selected)) : []),
      [selected]
    );
    const [name, setName] = React.useState<string>("");
    const [wasValidated, setWasValidated] = React.useState<boolean>(false);
    const [showStructures, setShowStructures] = React.useState<
      | {
          type: "selection";
        }
      | {
          type: "cluster";
          cluster: string;
        }
      | null
    >(null);
    const [rankingSelection, setRankingSelection] = React.useState<IParticleSelection>(null);
    const [heatmapData, setHeatmapData] = React.useState<
      { data: Partial<Plotly.PlotData>; layout: Partial<Plotly.Layout> }[] | null
    >(null);
    const [heatmapLoading, setHeatmapLoading] = React.useState<boolean>(false);
    const [heatmapHover, setHeatmapHover] = React.useState<string[] | null>(null);
    const [heatmapHoverIndices, setHeatmapHoverIndices] = React.useState<number[] | null>(null);

    const clusterNames = Object.keys(clusters);
    const clusterNameTaken = clusterNames.includes(name);

    const structuresToShow = React.useMemo<ICollection[] | null>(() => {
      if (showStructures) {
        const filteredParticles =
          showStructures.type === "selection"
            ? selectedParticles
            : showStructures.type === "cluster"
            ? clusters[showStructures.cluster]
            : null;
        return filteredParticles
          ? Object.entries(groupByLodash(filteredParticles, "collection")).map(([name, data]) => ({ name, data }))
          : null;
      }
      return null;
    }, [showStructures, clusters, selectedParticles]);

    return (
      <>
        <Modal
          show={Boolean(showStructures && structuresToShow)}
          onHide={() => setShowStructures(null)}
          size="xl"
          dialogClassName="modal-full-width"
        >
          {showStructures && structuresToShow ? (
            <>
              <Modal.Header closeButton>
                <Modal.Title>View {structuresToShow.reduce((a, b) => a + b.data.length, 0)} entries</Modal.Title>
              </Modal.Header>
              <Modal.Body
                style={{
                  display: "flex",
                  flexDirection: "column",
                }}
              >
                <StructureCardGrid collections={structuresToShow} setSelected={(s) => setRankingSelection(s)} />
              </Modal.Body>
              <Modal.Footer>
                <Button
                  variant="primary"
                  onClick={() => {
                    setShowStructures(null);
                    setSelected(rankingSelection);
                  }}
                  disabled={!rankingSelection}
                >
                  Save as new selection
                </Button>
                <Button variant="secondary" onClick={() => setShowStructures(null)}>
                  Close
                </Button>
              </Modal.Footer>
            </>
          ) : null}
        </Modal>
        <Modal
          show={Boolean(heatmapData)}
          onHide={() => setHeatmapData(null)}
          size="xl"
          dialogClassName="modal-full-width"
        >
          {heatmapData ? (
            <>
              <Modal.Header closeButton>
                <Modal.Title>View correlation</Modal.Title>
              </Modal.Header>
              <Modal.Body
                style={{
                  display: "flex",
                  flexDirection: "column",
                }}
              >
                <div style={{ height: 210, visibility: heatmapHover ? "visible" : "hidden", margin: "0px auto" }}>
                  {heatmapHover ? (
                    <div className="d-flex flex-row align-items-center" style={{ gap: 50 }}>
                      {heatmapHover.map((structure) => (
                        <StructureImage
                          structure={structure}
                          style={{
                            width: 210,
                          }}
                        />
                      ))}
                      {heatmapHover.length > 1 ? (
                        <StructureImage
                          structure={heatmapHover}
                          style={{
                            width: 210,
                          }}
                        />
                      ) : null}
                    </div>
                  ) : null}
                </div>
                <div className="d-flex flex-row">
                  {heatmapData.map(({ data, layout }, i) => (
                    <PlotComponent
                      style={{ height: "100%", flex: 1 }}
                      onHover={(e) => {
                        if (e.points[0]) {
                          setHeatmapHover(Array.from(new Set([e.points[0].x as string, e.points[0].y as string])));
                          setHeatmapHoverIndices((e.points[0].pointIndex as any) as number[]);
                        }
                      }}
                      onUnhover={() => {
                        setHeatmapHover(null);
                        setHeatmapHoverIndices(null);
                      }}
                      data={
                        heatmapHoverIndices
                          ? [
                              { ...data, opacity: 0.5 },
                              {
                                ...data,
                                showlegend: false,
                                showscale: false,
                                z: (() => {
                                  try {
                                    const originalZ = data.z! as number[][];
                                    // Copy the rows and cols only of the selected range
                                    if (
                                      heatmapHoverIndices &&
                                      heatmapHoverIndices.length === 2 &&
                                      originalZ.length > 0
                                    ) {
                                      const emptyZ = Array(originalZ.length)
                                        .fill(null)
                                        .map(() => Array(originalZ[0].length).fill(null));

                                      const x = heatmapHoverIndices[0],
                                        y = heatmapHoverIndices[1];
                                      emptyZ[x] = originalZ[x];
                                      emptyZ.forEach((row, i) => {
                                        row[y] = originalZ[i][y];
                                      });
                                      return emptyZ;
                                    }
                                  } catch {
                                    return data.z;
                                  }
                                })(),
                                opacity: 1,
                              },
                            ]
                          : [data]
                      }
                      layout={{
                        dragmode: "lasso",
                        hovermode: "closest",
                        autosize: true,
                        legend: {
                          x: 1,
                          xanchor: "right",
                          y: 1,
                        },
                        scene: {
                          aspectmode: "data",
                        },
                        margin: {
                          l: 250,
                          r: 0,
                          b: 0,
                          t: 25,
                          pad: 4,
                        },
                        xaxis: {
                          // automargin: true,
                        },
                        yaxis: {
                          scaleanchor: "x",
                          visible: i === 0,
                          autorange: "reversed",
                          // automargin: true,
                        },
                        ...layout,
                      }}
                      config={{...PLOTLY_CONFIG, displayModeBar: false}}
                    />
                  ))}
                </div>
              </Modal.Body>
              <Modal.Footer>
                <Button variant="secondary" onClick={() => setHeatmapData(null)}>
                  Close
                </Button>
              </Modal.Footer>
            </>
          ) : null}
        </Modal>
        {selected ? (
          <>
            <p className="lead">Selection</p>
            <p>
              {selectedParticles.length} Selected{" "}
              <div className="btn-group btn-group-sm ml-2 mr-2 float-right" role="group">
                <button
                  type="button"
                  className="btn btn-sm btn-light"
                  title="Show in table view"
                  onClick={() => setShowStructures({ type: "selection" })}
                >
                  <i className="fas fa-fw fa-table" />
                </button>
                <button
                  type="button"
                  className="btn btn-sm btn-light"
                  disabled={heatmapLoading || selectedParticles.length > 50}
                  title="Compute distances"
                  onClick={async () => {
                    setHeatmapLoading(true);
                    const selectedWithEmbeddings = await embedStructures({
                      structures: selectedParticles.map((p) => p.structure),
                      include_embedding: true,
                    });

                    if (selectedWithEmbeddings.length > 0) {
                      const cdddSimilarities = selectedWithEmbeddings.map((a) =>
                        selectedWithEmbeddings.map((b) =>
                          a.embedding!.reduce((acc, x, i) => {
                            const y = b.embedding![i];
                            return acc + Math.pow(x - y, 2);
                          }, 0)
                        )
                      );

                      const selectedStructures = selectedWithEmbeddings.map((p) => p.structure);
                      const tanimotoSimilaritiesRaw = await Promise.all(
                        selectedStructures.map((structure) =>
                          getTanimotoSimilarity(selectedStructures, structure, "morgan").then(
                            ({ tanimoto }) => tanimoto
                          )
                        )
                      );

                      const tanimotoSimilarities = selectedStructures.map((x, i) =>
                        selectedStructures.map((y) => tanimotoSimilaritiesRaw[i][y] || 0)
                      );

                      setHeatmapData(
                        [
                          {
                            data: { z: cdddSimilarities, reversescale: true, zmin: 0 },
                            layout: { title: "CDDD Distance" },
                          },
                          {
                            data: { z: tanimotoSimilarities, zmin: 0, zmax: 1 },
                            layout: { title: "Morgan Similarities" },
                          },
                        ].map(({ data, layout }) => ({
                          data: {
                            x: selectedStructures,
                            y: selectedStructures,
                            hoverinfo: "z",
                            colorscale: "YlGnBu",
                            type: "heatmap",
                            ...data,
                          },
                          layout,
                        }))
                      );
                    }
                    setHeatmapLoading(false);
                  }}
                >
                  <i className={`fas fa-fw ${heatmapLoading ? "fa-circle-notch fa-spin" : "fa-chess-board"}`} />
                </button>
              </div>
            </p>
          </>
        ) : null}
        {clusterNames.length > 0 ? (
          <>
            <p className="lead">Clusters</p>
            {clusterNames.map((cluster) => (
              <div key={cluster} className="clearfix mb-1">
                {cluster} ({clusters[cluster].length})
                <div className="float-right">
                  <div className="btn-group btn-group-sm ml-2 mr-2" role="group">
                    <button
                      type="button"
                      className="btn btn-light"
                      title="Show cluster in table view"
                      onClick={() => {
                        setShowStructures({ type: "cluster", cluster });
                      }}
                    >
                      <i className="fas fa-fw fa-table" />
                    </button>
                    <button
                      type="button"
                      className="btn btn-danger"
                      title="Delete cluster"
                      onClick={() => {
                        setClusters(cluster, null);
                      }}
                    >
                      <i className="fas fa-fw fa-times" />
                    </button>
                  </div>
                </div>
              </div>
            ))}
          </>
        ) : (
          <p className="lead">No clusters created</p>
        )}
        {selectedParticles.length > 0 ? (
          <>
            <p className="lead">Create new cluster</p>
            {selectedParticles.length === 1 ? (
              <StructureCard structure={selectedParticles[0]} showProperties={false} />
            ) : (
              <StructureImage style={{ maxWidth: "100%" }} structure={selectedParticles.map((s) => s.structure)} />
            )}
            <form
              className={`${wasValidated ? "was-validated" : ""}`}
              onSubmit={(e) => {
                e.preventDefault();
                e.stopPropagation();

                if (clusterNameTaken) {
                  setWasValidated(true);
                  return;
                }
                setClusters(name, selectedParticles);
                setWasValidated(false);
                setName("");
              }}
            >
              <div className="form-group">
                <label htmlFor="clusterNameInput">Cluster name</label>
                <input
                  type="text"
                  className={`form-control form-control-sm ${clusterNameTaken ? "is-invalid" : ""}`}
                  id="clusterNameInput"
                  value={name}
                  required={true}
                  onChange={(e) => setName(e.currentTarget.value)}
                />
                {clusterNameTaken ? (
                  <div className="invalid-feedback">Cluster with this name already exists.</div>
                ) : null}
              </div>

              <div className="text-right">
                <button type="submit" className="btn btn-primary" disabled={clusterNameTaken}>
                  Create cluster
                </button>
              </div>
            </form>
          </>
        ) : (
          <p className="lead">Select to create a new cluster</p>
        )}
      </>
    );
  }
);
