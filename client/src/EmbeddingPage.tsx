import * as React from "react";
import { LoadingPage } from "./components/LoadingPage";
import {
  DEFAULT_CHEMBL_COLLECTION,
  EActiveTabs,
  ICollection,
  IEnabledProjection,
  IParticle,
  IParticleSelection,
  IPlotOptions,
  IRegistry,
} from "./interfaces";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { HorizontalCollapse } from "./components/HorizontalCollapse";
import { ClusterSidePanel } from "./ClusterSidePanel";
import { Tooltip } from "./utils/hooks";
import { GridItemOptions } from "./components/GridItemOptions";
import { Grid } from "./components/Grid";
import { ComputeEmbeddingsForm, MSOForm, TanimotoForm, SubstructureMatchingForm } from "./components/form";
import { getChemblUMAPEmbedding } from "./utils/api";
import { ScatterPlot } from "./components/ScatterPlot";
import { InterpolationForm } from "./components/form/InterpolationForm";
import { MatchedMolecularPairsForm } from "./components/form/MatchedMolecularPairsForm";
import { NeighborSamplingForm } from "./components/form/NeighborSamplingForm";
import { RecomputeEmbeddingsForm } from "./components/form/RecomputeEmbeddingsForm";
import CreatableSelect from "react-select/creatable";
import { ParallelCoordinatesPlot } from "./components/ParallelCoordinatesPlot";
import { StonedSelfiesForm } from "./components/form/StonedSelfiesForm";
import { StructureImage } from "./components/StructureImage";
import { ProjectionSettingsModal } from "./components/ProjectionSettingsModal";
import { ProjectionSettingsForm } from "./components/form/ProjectionSettingsForm";
import pickBy from "lodash.pickby";
// @ts-ignore Typings?
import Split from "react-split";
import ReactResizeDetector from "react-resize-detector";

const NR_OF_LATENT_SPACE_PARTICLES = 10;

export function EmbeddingPage({
  registry,
  collections,
  setCollections,
  interpolationStructures,
  setInterpolationStructures,
  setActiveTab,
}: {
  registry: IRegistry | null;
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  interpolationStructures: string[];
  setInterpolationStructures(structures: string[]): void;
  setActiveTab(tab: EActiveTabs): void;
}) {
  const allParticles = React.useMemo(
    () => collections.reduce<IParticle[]>((acc, cur) => [...acc, ...cur.data], []),
    [collections]
  );

  const [visibleNeighborhoodSamplings, setVisibleNeighborhoodSamplings] = React.useState<IParticle[] | null>(null);
  const [loading, setLoading] = React.useState<boolean>(false);
  const [clusterCollapsed, setClusterCollapsed] = React.useState<boolean>(false);
  const [optionsCollapsed, setOptionsCollapsed] = React.useState<boolean>(false);
  const [enabledProjections, setEnabledProjections] = React.useState<IEnabledProjection[]>([
    {
      label: "chembl_umap",
      value: "chembl_umap",
      projection: "chembl_umap",
      plotOptions: {},
    },
  ]);
  const [editProjection, setEditProjection] = React.useState<IEnabledProjection | undefined>(undefined);

  const [hover, setHover] = React.useState<IParticle | null>(null);
  const [filtered, setFiltered] = React.useState<IParticleSelection>(null);
  // Wrap the setSelected into a shallow list compare to avoid rerenders
  // eslint-disable-next-line react-hooks/exhaustive-deps
  // const selected = React.useMemo(() => mainParticles.filter((p) => p.selected), [mainParticles, _selectedChanged]);
  // const selected = React.useMemo<{ [key: string]: IParticle[] } | undefined>(() => undefined, []);
  const [selection, _setSelection] = React.useState<IParticleSelection>(null);
  const setSelection = React.useCallback(
    (s: IParticleSelection) => {
      collections.forEach((c) => c.data.forEach((p) => (p.selected = false)));
      if (s) {
        const _setParticleSelected = (particle: IParticle) => (particle.selected = true);
        if (Array.isArray(s)) {
          s.forEach(_setParticleSelected);
        } else {
          Object.values(s).forEach((particles) => particles.forEach(_setParticleSelected));
        }
      }
      _setSelection(s);
    },
    [collections]
  );

  // Plot options
  const [plotOptions, setPlotOptions] = React.useState<IPlotOptions>({
    colorBy: null,
    opacityBy: null,
    groupBy: null,
    connectBy: null,
    sizeBy: null,
  });

  const [customX, setCustomX] = React.useState<string | null>(null);
  const [customY, setCustomY] = React.useState<string | null>(null);
  const [customPlotSettings, setCustomPlotSettings] = React.useState<IEnabledProjection | undefined>(undefined);

  // Ranking options
  const [showSelectedOnly, setShowSelectedOnly] = React.useState<boolean>(false);

  const [availableProperties, availableProjections, availableEmbeddings, availableNearestNeighbors, availableClusters] =
    React.useMemo(() => {
      const keys: (keyof IParticle)[] = ["properties", "projection", "embedding", "nearest_neighbors", "clusters"];
      return keys.map((field) =>
        Array.from(
          allParticles
            .reduce<Set<string>>((acc, cur) => {
              Object.keys(cur[field] || {}).forEach((key) => acc.add(key));
              return acc;
            }, new Set())
            .keys()
        )
      );
    }, [allParticles]);

  const availableOpacityProperties = React.useMemo(() => [...availableProperties, "constant"], [availableProperties]);
  const availableConnectByProperties = React.useMemo(() => availableProperties, [availableProperties]);

  React.useEffect(() => {
    if (plotOptions.colorBy && !availableProperties.includes(plotOptions.colorBy)) {
      // setPlotOptions({ ...plotOptions, colorBy: null });
    }
    if (
      plotOptions.opacityBy &&
      typeof plotOptions.opacityBy !== "number" &&
      !availableOpacityProperties.includes(plotOptions.opacityBy)
    ) {
      setPlotOptions({ ...plotOptions, opacityBy: null });
    }
    if (plotOptions.groupBy && !availableProperties.includes(plotOptions.groupBy)) {
      setPlotOptions({ ...plotOptions, groupBy: null });
    }
    if (
      plotOptions.connectBy &&
      plotOptions.connectBy.some((option) => !availableConnectByProperties.includes(option))
    ) {
      // setPlotOptions({ ...plotOptions, connectBy: null });
    }
    if (plotOptions.sizeBy && !availableProperties.includes(plotOptions.sizeBy)) {
      setPlotOptions({ ...plotOptions, sizeBy: null });
    }
  }, [availableProperties, availableOpacityProperties, availableConnectByProperties, plotOptions]);

  React.useEffect(
    () => {
      if (false && !collections.find((c) => c.name === DEFAULT_CHEMBL_COLLECTION)) {
        getChemblUMAPEmbedding()
          .then((data) => setCollections([...collections, { data, name: DEFAULT_CHEMBL_COLLECTION }]))
          .catch((e) => console.error("Error loading chembl umap", e));
      }
    },
    [
      /* collections, setCollections */
    ]
  );

  const filteredCollections = React.useMemo<ICollection[]>(() => {
    return filtered ? Object.entries(filtered).map(([name, data]) => ({ name, data })) : collections;
  }, [filtered, collections]);

  const visibleCollections = React.useMemo<ICollection[]>(() => {
    return showSelectedOnly && selection
      ? Object.entries(selection).map(([name, data]) => ({ name, data }))
      : collections;
  }, [collections, showSelectedOnly, selection]);

  const latentSpaceCollections = React.useMemo<{ [key: string]: ICollection[] } | null>(() => {
    return true
      ? availableEmbeddings.reduce<{ [key: string]: any }>(
          (acc, cur) => ({
            ...acc,
            [cur]: filteredCollections.map((c) => {
              // Restrict to the first n particles
              const data = c.data
                .slice(0, NR_OF_LATENT_SPACE_PARTICLES)
                .filter((d) => d.embedding?.[cur])
                .map((d) => {
                  return d.embedding![cur].map((y, x) => ({
                    ...d,
                    properties: {
                      ...(d.properties || {}),
                      latentX: x,
                      latentY: y,
                    },
                  }));
                })
                .flat();
              return {
                ...c,
                data,
              };
            }),
          }),
          {}
        )
      : null;
  }, [filteredCollections, availableEmbeddings]);

  return (
    <>
      <Tooltip>
        {hover ? (
          // <StructureCard
          //   structure={hover}
          //   showProperties={false}
          //   style={{
          //     width: 180,
          //     backgroundColor: 'rgba(255, 255, 255, 0.2)'
          //   }}
          // />
          <StructureImage
            structure={hover.structure}
            style={{
              width: 180,
              backgroundColor: "rgba(255, 255, 255, 0.6)",
            }}
          />
        ) : null}
      </Tooltip>
      <ProjectionSettingsModal
        config={editProjection}
        availableNearestNeighbors={availableNearestNeighbors}
        availableClusters={availableClusters}
        availableProjections={availableProjections}
        availableProperties={availableProperties}
        availableOpacityProperties={availableOpacityProperties}
        availableConnectByProperties={availableConnectByProperties}
        onHide={() => setEditProjection(undefined)}
        onSave={(value) => {
          setEditProjection(undefined);
          // TODO: Hacky way to distinguish if the custom plot or one of the added plots was edited.
          if(enabledProjections.includes(editProjection!)) {
            setEnabledProjections((projections) => projections.map((p) => (p === editProjection ? value : p)));
          } else if(!customPlotSettings || customPlotSettings === editProjection) {
            setCustomPlotSettings(value);
          }
        }}
      />
      <HorizontalCollapse
        label="Options"
        position="left"
        size="col-3"
        collapsed={optionsCollapsed}
        setCollapsed={setOptionsCollapsed}
      >
        <ComputeEmbeddingsForm
          addCollection={(collection) => setCollections([...collections, collection])}
          loading={loading}
          setLoading={setLoading}
        />
        <MSOForm
          availableObjectives={registry?.objectives}
          addCollection={(collection) => setCollections([...collections, collection])}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
        />
        <InterpolationForm
          open={false}
          setCollection={(collection) => {
            setCollections([...collections.filter((c) => c.name !== collection.name), collection]);
          }}
          loading={loading}
          setLoading={setLoading}
          setStructures={setInterpolationStructures}
          structures={interpolationStructures}
        />
        <MatchedMolecularPairsForm
          addCollection={(collection) => setCollections([...collections, collection])}
          loading={loading}
          setLoading={setLoading}
          selection={selection}
        />
        <SubstructureMatchingForm
          collections={collections}
          setCollections={setCollections}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
        />
        <TanimotoForm
          collections={collections}
          setCollections={setCollections}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
        />
        <NeighborSamplingForm
          visible={visibleNeighborhoodSamplings}
          setVisible={setVisibleNeighborhoodSamplings}
          collections={collections}
          setCollections={setCollections}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
        />
        <RecomputeEmbeddingsForm
          collections={collections}
          setCollections={setCollections}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
        />
        <StonedSelfiesForm
          addCollection={(collection) => setCollections([...collections, collection])}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
        />
      </HorizontalCollapse>
      {/*<div className="mt-5" style={{ position: "relative", width: 0 }}>
        <div className="sticky-top" style={{ top: "20%", width: 100 }}>
          {hover ? (
            <StructureCard
              structure={hover}
              style={{
                left: 5,
                width: "100%",
              }}
            />
          ) : null}
        </div>
      </div>*/}
      <Split
        // minSize={250}
        gutterSize={5}
        sizes={[75, 25]}
        className={`split col-${optionsCollapsed ? "12" : "9"}`}
        style={{
          height: "100%",
          paddingLeft: optionsCollapsed ? 40 : 0,
          // marginRight: 40,
        }}
        onDragEnd={() => {
          window.dispatchEvent(new Event('resize'));
        }}
      >
        <div
          // className="col-md-10"
          style={{
            display: "flex",
            flexDirection: "column",
            height: "100%",
            overflow: "auto",
            paddingRight: 15,
          }}
        >
          <LoadingPage
            loading={loading && collections.length === 0}
            fallback="Please select structures for embedding"
            loadingText="Computing embedding of structures..."
          >
            {collections.length > 0 ? (
              <>
                <details>
                  <summary>
                    <h6 className="d-inline-block">
                      Plot Options <small className="text-muted">to customize the projection plots</small>
                    </h6>
                  </summary>
                  <div className="ms-3">
                    {[
                      {
                        name: "General",
                        plotOptions,
                        setter: (options: Partial<IPlotOptions>) => setPlotOptions({ ...plotOptions, ...options }),
                      },
                      ...collections.map((collection) => ({
                        name: collection.name,
                        plotOptions: collection.plotOptions || {},
                        setter: (options: Partial<IPlotOptions>) =>
                          setCollections(
                            collections.map((c) =>
                              c === collection
                                ? {
                                    ...c,
                                    plotOptions: {
                                      ...(c.plotOptions || {}),
                                      ...options,
                                    },
                                  }
                                : c
                            )
                          ),
                      })),
                    ].map(({ name, setter, plotOptions }, i) => (
                      <details open={i === 0}>
                        <summary>
                          {name} Options{" "}
                          <small className="text-muted">{Object.values(plotOptions).filter(Boolean).length}</small>
                        </summary>
                        <ProjectionSettingsForm
                          inline={true}
                          availableNearestNeighbors={availableNearestNeighbors}
                          availableClusters={availableClusters}
                          availableConnectByProperties={availableConnectByProperties}
                          availableProperties={availableProperties}
                          availableOpacityProperties={availableOpacityProperties}
                          setConfig={(newConfig) => setter(newConfig.plotOptions)}
                          config={{
                            label: "General",
                            value: "",
                            projection: "",
                            plotOptions,
                          }}
                        />
                      </details>
                    ))}
                  </div>
                </details>

                <details className="row" style={{ position: "relative" }} open={true}>
                  <summary>
                    <h6 className="d-inline-block">
                      Parallel Coordinates <small className="text-muted">to compare multiple numeric properties</small>
                    </h6>
                  </summary>
                  <ParallelCoordinatesPlot
                    collections={filteredCollections}
                    selection={selection}
                    setSelection={setSelection}
                  />
                </details>

                <details className="row" style={{ position: "relative" }} open={true}>
                  <summary>
                    <h6 className="d-inline-block">
                      Custom Plot <small className="text-muted">to compare two arbitrary properties</small>
                    </h6>
                  </summary>
                  <Grid>
                    <GridItemOptions
                      key="Custom"
                      gridOptions={{
                        w: 12,
                        h: 25,
                        y: 0,
                        x: 0,
                      }}
                      onSettings={() => setEditProjection(customPlotSettings || {
                        label: 'Custom',
                        value: 'Custom',
                        projection: null,
                        plotOptions: {}
                      })}
                    >
                      <ScatterPlot
                        title={customX && customY ? `${customX} vs. ${customY}` : "Please select two axis properties"}
                        collections={filteredCollections}
                        options={{ ...plotOptions, ...pickBy(customPlotSettings?.plotOptions || {}, (v) => v) }}
                        hover={hover}
                        setHover={setHover}
                        selected={selection}
                        setSelected={setSelection}
                        xAccessor={`properties.${customX}`}
                        yAccessor={`properties.${customY}`}
                      >
                        <div
                          style={{
                            width: 150,
                            position: "absolute",
                            left: -45,
                            top: "50%",
                            transform: "rotate(-90deg)",
                          }}
                        >
                          <select
                            className="form-control form-control-sm"
                            value={customY || ""}
                            onChange={(e) => setCustomY(e.currentTarget.value || null)}
                          >
                            <option value="">Select...</option>
                            {availableProperties.map((p) => (
                              <option key={p}>{p}</option>
                            ))}
                          </select>
                        </div>
                        <div
                          style={{
                            width: 150,
                            position: "absolute",
                            left: "40%",
                            bottom: 0,
                          }}
                        >
                          <select
                            className="form-control form-control-sm"
                            value={customX || ""}
                            onChange={(e) => setCustomX(e.currentTarget.value || null)}
                          >
                            <option value="">Select...</option>
                            {availableProperties.map((p) => (
                              <option key={p}>{p}</option>
                            ))}
                          </select>
                        </div>
                      </ScatterPlot>
                    </GridItemOptions>
                  </Grid>
                </details>

                <details className="row" style={{ position: "relative" }} open={true}>
                  <summary>
                    <h6 className="d-inline-block">
                      Latent Space Visualization <small className="text-muted"></small>
                    </h6>
                  </summary>
                  <Grid>
                    {Object.keys(latentSpaceCollections || []).map((embedding) => (
                      <GridItemOptions
                        key={embedding}
                        gridOptions={{
                          w: 12,
                          h: 25,
                          y: 0,
                          x: 0,
                        }}
                      >
                        <ScatterPlot
                          title={`Top ${NR_OF_LATENT_SPACE_PARTICLES} in ${embedding.toUpperCase()} Latent Space`}
                          collections={latentSpaceCollections![embedding]}
                          options={plotOptions}
                          hover={hover}
                          setHover={(hover) => {
                            setHover(
                              hover
                                ? collections
                                    .find((c) => c.name === hover.collection)
                                    ?.data.find((d) => d.index === hover.index) || null
                                : null
                            );
                          }}
                          selected={selection}
                          setSelected={(selected) => {
                            setSelection(selected);
                          }}
                          xAccessor={`properties[latentX]`}
                          yAccessor={`properties[latentY]`}
                        />
                      </GridItemOptions>
                    ))}
                  </Grid>
                </details>

                <details open={true}>
                  <summary>
                    <h6 className="d-inline-block">
                      Projections{" "}
                      <small className="text-muted">
                        {enabledProjections.length} out of {availableProjections.length}
                      </small>
                    </h6>
                  </summary>
                  <CreatableSelect<IEnabledProjection, true>
                    isMulti
                    name="projections"
                    value={enabledProjections}
                    onChange={(e) => {
                      setEnabledProjections(e.map((v) => v));
                    }}
                    onCreateOption={(label) => {
                      const newOption: IEnabledProjection = {
                        label: label,
                        value: label,
                        projection: "",
                        plotOptions: {},
                      };
                      setEnabledProjections([...enabledProjections, newOption]);
                      setEditProjection(newOption);
                    }}
                    options={availableProjections.map((p) => ({
                      label: p,
                      value: p,
                      projection: p,
                      plotOptions: {},
                    }))}
                    closeMenuOnSelect={false}
                  />
                  <div className="row" style={{ position: "relative" }}>
                    <Grid>
                      {
                        enabledProjections
                          // .filter((p) => availableProjections.includes(p.projection))
                          .map((config, i) => (
                            <GridItemOptions
                              key={config.value}
                              onClose={() => {
                                setEnabledProjections(enabledProjections.filter((p) => p !== config));
                              }}
                              onSettings={() => setEditProjection(config)}
                              gridOptions={{
                                w: 6,
                                h: 25,
                                y: 0 + Math.floor(i / 2) * 25,
                                x: i % 2 === 0 ? 0 : 6,
                              }}
                              renderInfo={() => (
                                <>
                                  {filteredCollections.map((c) => (
                                    <div className="d-flex flex-column">
                                      <strong className="text-nowrap text-truncate">{c.name}</strong>
                                      {Object.entries(c.projections?.[config.projection!] || {})
                                        .filter(([key, value]) => typeof value === "number")
                                        .map(([key, value]) => (
                                          <div className="text-nowrap d-flex">
                                            <span className="text-truncate">{key}</span>
                                            <span className="flex-fill">
                                              : {(value as number).toFixed(3) || "Not available"}
                                            </span>
                                          </div>
                                        ))}
                                    </div>
                                  ))}
                                </>
                              )}
                            >
                              <ScatterPlot
                                title={config.value}
                                collections={filteredCollections}
                                // Pick only "valid" values, as otherwise null would override the previous value
                                options={{ ...plotOptions, ...pickBy(config.plotOptions, (v) => v) }}
                                hover={hover}
                                setHover={setHover}
                                selected={selection}
                                setSelected={setSelection}
                                xAccessor={`projection[${config.projection}][0]`}
                                yAccessor={`projection[${config.projection}][1]`}
                              />
                            </GridItemOptions>
                          )) as any
                      }
                    </Grid>
                  </div>
                </details>

                <StructureCardGrid
                  collections={visibleCollections}
                  selected={selection}
                  setSelected={setSelection}
                  setFiltered={setFiltered}
                  tableClass="main-ranking"
                  structureCardProps={(structure) => ({
                    className: structure === hover ? "border-primary" : structure.selected ? "border-secondary" : "",
                  })}
                  renderTopForm={
                    <>
                      <div className="form-check form-switch me-sm-2">
                        <input
                          type="checkbox"
                          className="form-check-input"
                          id="customSwitch1"
                          checked={showSelectedOnly}
                          onChange={(e) => setShowSelectedOnly(e.currentTarget.checked)}
                        />
                        <label className="form-check-label" htmlFor="customSwitch1">
                          Show selected only
                        </label>
                      </div>
                      <button
                        className="btn btn-sm btn-primary"
                        disabled={!selection}
                        onClick={() => {
                          setInterpolationStructures(
                            Object.values(selection!)
                              .flat()
                              .map(({ structure }) => structure)
                          );
                          setActiveTab(EActiveTabs.INTERPOLATION);
                        }}
                      >
                        Use for interpolation
                      </button>
                    </>
                  }
                />
              </>
            ) : null}
          </LoadingPage>
        </div>
        <div>
          {/* <HorizontalCollapse
        label="Collections"
        position="right"
        size="col-md-4"
        collapsed={clusterCollapsed}
        setCollapsed={setClusterCollapsed}
      > */}
          <ReactResizeDetector handleWidth>
            {({ width, height }) => (
              <div
                style={{
                  height: "100%",
                  overflow: "auto",
                  paddingRight: 15,
                }}
              >
                {width! > 200 ? (
                  <ClusterSidePanel
                    collections={collections}
                    setCollections={setCollections}
                    selection={selection}
                    setSelection={setSelection}
                    setVisibleNeighborhoodSamplings={setVisibleNeighborhoodSamplings}
                    hover={hover}
                    setHover={setHover}
                  />
                ) : (
                  <span style={{ transform: "rotate(90deg)", position: "absolute", top: "50%", whiteSpace: "nowrap" }}>
                    Additional Options
                  </span>
                )}
              </div>
            )}
          </ReactResizeDetector>
          {/* </HorizontalCollapse> */}
        </div>
      </Split>
    </>
  );
}
