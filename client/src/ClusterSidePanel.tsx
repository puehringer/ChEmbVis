import * as React from "react";
import { Modal, Button } from "react-bootstrap";
import { StructureCard } from "./components/StructureCard";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { StructureImage } from "./components/StructureImage";
import { ICollection, IParticle, IParticleSelection, IServerCollection } from "./interfaces";
import groupByLodash from "lodash.groupby";
import { embedStructures, getTanimotoSimilarity } from "./utils/api";
import { PlotComponent, PLOTLY_CONFIG } from "./components/PlotComponent";
import { useNameInput } from "./utils/hooks";
import { FormWrapper } from "./components/form";
import { ButtonWithUpload } from "./components/ButtonWithUpload";
import cloneDeepWith from "lodash.clonedeepwith";
import { ParallelCoordinatesPlot } from "./components/ParallelCoordinatesPlot";
import { downloadCSVFile, normalizeArray } from "./utils";
import { Tabs, Tab } from "react-bootstrap";
import { extent } from "d3-array";
import { LocalNeighborhoodPlot } from "./components/LocalNeighborhoodPlot";
import Select from "react-select/src/Select";
import { GroupFlowSankeyPlot } from "./components/GroupFlowSankeyPlot";

export const ClusterSidePanel = React.memo(
  ({
    collections,
    setCollections,
    selection,
    setSelection,
    hover,
    setHover,
    setVisibleNeighborhoodSamplings,
  }: {
    collections: ICollection[];
    setCollections(collections: ICollection[]): void;
    selection: IParticleSelection;
    setSelection(selected: IParticleSelection): void;
    hover: IParticle | null;
    setHover(hover: IParticle | null): void;
    setVisibleNeighborhoodSamplings(collection: IParticle[] | null): void;
  }) => {
    const selectedParticles = React.useMemo<IParticle[]>(
      () => (selection ? ([] as IParticle[]).concat(...Object.values(selection)) : []),
      [selection]
    );

    const hoverOrSelected = hover || (selectedParticles.length === 1 ? selectedParticles[0] : null);
    const selectedOrHover = (selectedParticles.length === 1 ? selectedParticles[0] : null) || hover;

    const [lastHoverOrSelected, setHoverOrLastSelected] = React.useState<IParticle | null>(hoverOrSelected);

    React.useEffect(() => {
      if (hoverOrSelected) {
        setHoverOrLastSelected(hoverOrSelected);
      }
    }, [hoverOrSelected]);

    const selectionCollection =
      selection && Object.entries(selection).length === 1
        ? collections.find((c) => c.name === Object.keys(selection)[0])
        : null;
    const [visibleProjectionModels, setVisibleProjectionModels] = React.useState<
      IServerCollection["projections"] | null
    >(null);

    const [name, setName, nameInput] = useNameInput("newCollectionNameInput", "");
    const [recomputeEmbeddings, setRecomputeEmbeddings] = React.useState<boolean>(false);
    const [knnEmbeddingKey, setKnnEmbeddingKey] = React.useState<string>("");
    const [removeFromOriginal, setRemoveFromOriginal] = React.useState<boolean>(false);
    const [includeNeighbors, setIncludeNeighbors] = React.useState<{ [key: string]: number } | null>(null);
    const [newCollectionLoading, setNewCollectionLoading] = React.useState<boolean>(false);
    const [showStructures, setShowStructures] = React.useState<
      | {
          type: "selection";
        }
      | {
          type: "collection";
          collection: ICollection;
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

    const [parallelCoordinatesCollections, setParallelCoordinatesCollections] = React.useState<ICollection[] | null>(
      null
    );

    const structuresToShow = React.useMemo<ICollection[] | null>(() => {
      if (showStructures) {
        if (showStructures.type === "collection") {
          return [showStructures.collection];
        }
        const filteredParticles = showStructures.type === "selection" ? selectedParticles : null;
        return filteredParticles
          ? Object.entries(groupByLodash(filteredParticles, "collection")).map(([name, data]) => ({ name, data }))
          : null;
      }
      return null;
    }, [showStructures, selectedParticles]);

    return (
      <>
        <Modal show={Boolean(visibleProjectionModels)} onHide={() => setVisibleProjectionModels(null)} size="xl">
          <Modal.Header closeButton>
            <Modal.Title>Available projection models</Modal.Title>
          </Modal.Header>
          <Modal.Body>
            {visibleProjectionModels ? (
              <ul>
                {Object.entries(visibleProjectionModels)
                  .filter(([name, value]) => value.model)
                  .map(([name, value]) => (
                    <li key={name} className="text-truncate">
                      {name}: <span title={value.model || ""}>{value.model}</span>
                    </li>
                  ))}
              </ul>
            ) : null}
          </Modal.Body>
        </Modal>
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
                    setSelection(rankingSelection);
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
                          setHeatmapHoverIndices(e.points[0].pointIndex as any as number[]);
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
                      config={{ ...PLOTLY_CONFIG, displayModeBar: false }}
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
        <Modal
          show={Boolean(parallelCoordinatesCollections)}
          onHide={() => setParallelCoordinatesCollections(null)}
          size="xl"
          dialogClassName="modal-full-width"
        >
          {parallelCoordinatesCollections ? (
            <>
              <Modal.Header closeButton>
                <Modal.Title>
                  View parallel coordinates of {parallelCoordinatesCollections.map((c) => c.name).join(", ")}
                </Modal.Title>
              </Modal.Header>
              <Modal.Body>
                <ParallelCoordinatesPlot
                  collections={parallelCoordinatesCollections}
                  selection={selection}
                  setSelection={setSelection}
                />
              </Modal.Body>
              <Modal.Footer>
                <Button variant="secondary" onClick={() => setParallelCoordinatesCollections(null)}>
                  Close
                </Button>
              </Modal.Footer>
            </>
          ) : null}
        </Modal>
        <Tabs defaultActiveKey="collections" className="mt-3 mb-3" unmountOnExit={true}>
          <Tab eventKey="collections" title="Collections">
            {collections.length > 0 ? (
              <>
                {collections.map((c) => {
                  return (
                    <div key={c.name} className="d-flex align-items-center mb-1">
                      <span className="text-truncate me-auto" title={c.name}>
                        {c.name} <small>{c.data.length}</small>
                      </span>
                      <div className="btn-group btn-group-sm ms-2 me-2" role="group">
                        <button
                          type="button"
                          className="btn btn-light"
                          title="Show parallel coordinates"
                          onClick={() => {
                            setParallelCoordinatesCollections([c]);
                          }}
                        >
                          <i className="fas fa-fw fa-arrows-alt-h" />
                        </button>
                        {Object.entries(c.projections || {}).some(([key, value]) => value.model) ? (
                          <button
                            type="button"
                            className="btn btn-light"
                            title="Show attached projection models"
                            onClick={() => {
                              setVisibleProjectionModels(c.projections);
                            }}
                          >
                            <i className="fas fa-fw fa-server" />
                          </button>
                        ) : null}
                        {c.type === "neighborhoodSampling" ? (
                          <button
                            type="button"
                            className="btn btn-light"
                            title="Show neighborhood grid"
                            onClick={() => {
                              setVisibleNeighborhoodSamplings(c.data);
                            }}
                          >
                            <i className="fas fa-fw fa-th" />
                          </button>
                        ) : null}
                        <button
                          type="button"
                          className="btn btn-light"
                          title="Download SMILES as CSV"
                          onClick={() => {
                            downloadCSVFile(c.data.map((d) => d.structure).join("\n"), "export");
                          }}
                        >
                          <i className="fas fa-fw fa-cloud-download-alt" />
                        </button>
                        <button
                          type="button"
                          className="btn btn-light"
                          title="Show collection in table view"
                          onClick={() => {
                            setShowStructures({
                              type: "collection",
                              collection: c,
                            });
                          }}
                        >
                          <i className="fas fa-fw fa-table" />
                        </button>
                        <button
                          type="button"
                          className="btn btn-light"
                          title={c.hidden ? "Enable collection" : "Disable collection"}
                          onClick={() => {
                            setCollections(
                              collections.map((collection) =>
                                c === collection ? { ...c, hidden: !c.hidden } : collection
                              )
                            );
                          }}
                        >
                          <i className={`fas fa-fw ${c.hidden ? "fa-eye-slash" : "fa-eye"}`} />
                        </button>
                        <button
                          type="button"
                          className="btn btn-danger"
                          title="Delete collection"
                          onClick={() => {
                            setCollections(collections.filter((collection) => c !== collection));
                          }}
                        >
                          <i className="fas fa-fw fa-times" />
                        </button>
                      </div>
                    </div>
                  );
                })}
              </>
            ) : (
              <p className="lead">No collections found</p>
            )}
          </Tab>
          <Tab eventKey="selection" title="Selection">
            {selection ? (
              <>
                <p>
                  {selectedParticles.length} Selected{" "}
                  <div className="btn-group btn-group-sm ms-2 me-2 float-end" role="group">
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
                        const selectedWithEmbeddings = (
                          await embedStructures({
                            structures: selectedParticles.map((p) => p.structure),
                            include_embedding: true,
                          })
                        ).data;

                        if (selectedWithEmbeddings.length > 0) {
                          const cdddSimilarities = selectedWithEmbeddings.map((a) =>
                            selectedWithEmbeddings.map((b) =>
                              a.embedding!["cddd"].reduce((acc, x, i) => {
                                const y = b.embedding!["cddd"][i];
                                return acc + Math.pow(x - y, 2);
                              }, 0)
                            )
                          );

                          const selectedStructures = selectedWithEmbeddings.map((p) => p.structure);
                          const tanimotoSimilaritiesRaw = await Promise.all(
                            selectedStructures.map((structure) =>
                              getTanimotoSimilarity(selectedStructures, structure, "ecfp4").then(
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
            {selectedParticles.length > 0 ? (
              <FormWrapper
                title="Create new collection from selection"
                open={true}
                loading={newCollectionLoading}
                setLoading={setNewCollectionLoading}
                onSubmit={async () => {
                  if (selection) {
                    let data = Object.values(selection).flat();

                    if (includeNeighbors && data.length === 1) {
                      const neighbors = Object.entries(includeNeighbors)
                        .filter(([key, value]) => value)
                        .map(([key, value]) => data[0].nearest_neighbors?.[key]?.knn_particles?.slice(0, value)!)
                        .flat()
                        .filter(Boolean);
                      data.push(...neighbors);
                      data = Array.from(new Set(data));
                    }

                    const newCollection = recomputeEmbeddings
                      ? {
                          name,
                          ...(await embedStructures({
                            structures: data.map((d) => d.structure),
                            include_embedding: true,
                          })),
                        }
                      : {
                          name,
                          data: data.map((d) => ({
                            ...cloneDeepWith(d, (value, key, object, stack) =>
                              key === "knn_particles" ? undefined : value
                            ),
                            nearest_neighbors: undefined,
                          })),
                          // Add the projection model if only a single collection was selected
                          projections: selectionCollection?.projections,
                        };
                    const oldCollections = removeFromOriginal
                      ? collections.map((c) =>
                          selection[c.name]
                            ? {
                                ...c,
                                data: c.data.filter((d) => !selection[c.name].includes(d)),
                              }
                            : c
                        )
                      : collections;
                    setCollections([...oldCollections, newCollection]);
                    setSelection(null);
                    setName("");
                  }
                }}
              >
                <div className="d-flex justify-content-center">
                  {selectedParticles.length === 1 ? (
                    <StructureCard structure={selectedParticles[0]} showProperties={false} />
                  ) : (
                    <StructureImage
                      style={{ maxWidth: "100%" }}
                      structure={selectedParticles.map((s) => s.structure)}
                    />
                  )}
                </div>
                {nameInput}
                <div className="form-check form-switch me-sm-2">
                  <input
                    type="checkbox"
                    className="form-check-input"
                    id="collectionsRecomputeEmbeddingsInput"
                    checked={recomputeEmbeddings}
                    onChange={(e) => setRecomputeEmbeddings(e.currentTarget.checked)}
                  />
                  <label className="form-check-label" htmlFor="collectionsRecomputeEmbeddingsInput">
                    Recompute embeddings, projections, ...
                  </label>
                </div>
                <div className="form-check form-switch me-sm-2">
                  <input
                    type="checkbox"
                    className="form-check-input"
                    id="collectionsRemoveOldInput"
                    checked={removeFromOriginal}
                    onChange={(e) => setRemoveFromOriginal(e.currentTarget.checked)}
                  />
                  <label className="form-check-label" htmlFor="collectionsRemoveOldInput">
                    Remove from original collection
                  </label>
                </div>
                {selectedParticles.length === 1 && selectedParticles[0].nearest_neighbors ? (
                  <div>
                    <div className="form-check form-switch me-sm-2">
                      <input
                        type="checkbox"
                        className="form-check-input"
                        id="collectionsEnableKNN"
                        checked={Boolean(includeNeighbors)}
                        onChange={(e) => setIncludeNeighbors(e.currentTarget.checked ? {} : null)}
                      />
                      <label className="form-check-label" htmlFor="collectionsEnableKNN">
                        Enable KNN
                      </label>
                    </div>
                    {includeNeighbors &&
                      Object.entries(selectedParticles[0].nearest_neighbors || {}).map(([key, value]) => (
                        <div className="mb-3">
                          <label htmlFor="neighborhoodSamplesInput" className="form-label">
                            {key} {includeNeighbors?.[key]}
                          </label>
                          <input
                            type="range"
                            className="form-range"
                            id="neighborhoodSamplesInput"
                            value={includeNeighbors?.[key] || 0}
                            onChange={(e) =>
                              setIncludeNeighbors({ ...(includeNeighbors || {}), [key]: e.currentTarget.valueAsNumber })
                            }
                            min={0}
                            max={value.knn_ind.length}
                            step={1}
                          />
                        </div>
                      ))}
                  </div>
                ) : null}
                <div className="text-end">
                  <ButtonWithUpload loading={newCollectionLoading} disabled={!selection || !name} text="Create" />
                </div>
              </FormWrapper>
            ) : (
              <p className="lead">Select to create a new collection</p>
            )}
          </Tab>
          <Tab eventKey="knn" title="KNN">
            {true || hoverOrSelected ? (
              <>
                <LocalNeighborhoodPlot selected={lastHoverOrSelected} collection={collections.find((c) => c.name === hoverOrSelected?.collection)} setHover={setHover} hover={hover} />
                {/*<select
                  value={knnEmbeddingKey}
                  onChange={(e) => setKnnEmbeddingKey(e.currentTarget.value)}
                  className="custom-select custom-select-sm mb-1"
                >
                  <option value="">Choose embedding...</option>
                  {Object.entries(hoverOrSelected.nearest_neighbors || {}).map(([key, value]) => (
                    <option value={key}>{key}</option>
                  ))}
                </select>
                {knnEmbeddingKey ? (<>
                  <div
                    style={{
                      display: "grid",
                      gridTemplateColumns: "1fr 1fr 1fr 1fr",
                      gridColumnGap: 3,
                      gridRowGap: 3,
                    }}
                  >
                    {hoverOrSelected.nearest_neighbors?.[knnEmbeddingKey]?.knn_particles.map((p) => (
                      <div className="d-flex overflow-auto">
                        <StructureImage style={{ width: "100%" }} structure={p.structure} />
                      </div>
                    ))}
                  </div>
                </>) : null} */}
              </>
            ) : null}
          </Tab>
          <Tab eventKey="sankey" title="Group Flow">
            <GroupFlowSankeyPlot selection={selection} collections={collections} />
          </Tab>
          <Tab eventKey="details" title="Details">
            {hoverOrSelected ? (
              <StructureCard structure={hoverOrSelected} />
            ) : (
              <p className="lead">Hover over a structure to see details</p>
            )}
          </Tab>
        </Tabs>
      </>
    );
  }
);
