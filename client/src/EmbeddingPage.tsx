import * as React from "react";
import { LoadingPage } from "./components/LoadingPage";
import { PlotSelect } from "./components/PlotSelect";
import { StructureCard } from "./components/StructureCard";
import { EActiveTabs, ICollection, IParticle, IParticleSelection, IPlotOptions } from "./interfaces";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { HorizontalCollapse } from "./components/HorizontalCollapse";
import { ClusterSidePanel, getClustersFromParticle } from "./ClusterSidePanel";
import { Tooltip } from "./utils/hooks";
import { GridItemOptions } from "./components/GridItemOptions";
import { Grid } from "./components/Grid";
import { ComputeEmbeddingsForm, MSOForm, TanimotoForm, SubstructureMatchingForm } from "./components/form";
import { getChemblUMAPEmbedding } from "./utils/api";
import { ScatterPlot } from "./components/ScatterPlot";

export function EmbeddingPage({
  collections,
  setCollections,
  setInterpolationStructures,
  setActiveTab,
}: {
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  setInterpolationStructures(structures: string[]): void;
  setActiveTab(tab: EActiveTabs): void;
}) {
  const allParticles = React.useMemo(() => collections.reduce<IParticle[]>((acc, cur) => [...acc, ...cur.data], []), [
    collections,
  ]);

  const [loading, setLoading] = React.useState<boolean>(false);
  const [clusterCollapsed, setClusterCollapsed] = React.useState<boolean>(true);
  const [optionsCollapsed, setOptionsCollapsed] = React.useState<boolean>(false);
  const [enabledProjections, setEnabledProjections] = React.useState<string[]>(["chembl_umap"]);

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

  // Wrap the setClustersChanged into a shallow list compare to avoid rerenders
  const [_clustersChanged, _setClustersChanged] = React.useState<any>({});
  // eslint-disable-next-line react-hooks/exhaustive-deps
  // const clusters = React.useMemo(() => particles.filter((p) => p.selected), [
  //   particles,
  //   _selectedChanged,
  // ]);
  const clusters = React.useMemo<{ [key: string]: IParticle[] }>(() => {
    void _clustersChanged; // Line only exists to avoid eslint-disable-next-line react-hooks/exhaustive-deps
    const c: { [key: string]: IParticle[] } = {};
    allParticles.forEach((p) => {
      getClustersFromParticle(p).forEach((cluster) => {
        if (!c[cluster]) {
          c[cluster] = [];
        }
        c[cluster].push(p);
      });
    });
    // Sort alphabetically
    return Object.keys(c)
      .sort()
      .reduce((acc, cur) => ({ ...acc, [cur]: c[cur] }), {});
  }, [allParticles, _clustersChanged]);

  const setClusters = React.useCallback(
    (name: string, cluster: IParticle[] | null) => {
      if (!cluster) {
        allParticles.forEach((p) => {
          const c = getClustersFromParticle(p);
          if (c.includes(name)) {
            p.properties!.clusters = c
              .filter((clusterName) => clusterName !== name)
              .sort()
              .join(";");
          }
        });
      } else {
        cluster.forEach((p) => {
          p.properties = {
            ...(p.properties || {}),
            clusters: `${p.properties?.clusters ? `${p.properties.clusters};` : ""}${name}`,
          };
        });
      }
      _setClustersChanged({});
      setCollections([...collections]);
    },
    [allParticles, setCollections, collections]
  );

  // Plot options
  const [plotOptions, setPlotOptions] = React.useState<IPlotOptions>({
    colorBy: null,
    opacityBy: null,
    groupBy: null,
    connectByValues: null,
    sizeBy: null,
  });
  const [customX, setCustomX] = React.useState<string | null>(null);
  const [customY, setCustomY] = React.useState<string | null>(null);

  // Ranking options
  const [showSelectedOnly, setShowSelectedOnly] = React.useState<boolean>(false);

  const [availableProperties, availableProjections] = React.useMemo(() => {
    void clusters; // Line only exists to avoid eslint-disable-next-line react-hooks/exhaustive-deps
    const keys: (keyof IParticle)[] = ["properties", "projection"];
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
  }, [allParticles, clusters]);
  const availableOpacityProperties = React.useMemo(() => [...availableProperties, "constant"], [availableProperties]);

  React.useEffect(() => {
    if (plotOptions.colorBy && !availableProperties.includes(plotOptions.colorBy)) {
      setPlotOptions({ ...plotOptions, colorBy: null });
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
      plotOptions.connectByValues &&
      plotOptions.connectByValues.some((option) => !availableProperties.includes(option))
    ) {
      setPlotOptions({ ...plotOptions, connectByValues: null });
    }
    if (plotOptions.sizeBy && !availableProperties.includes(plotOptions.sizeBy)) {
      setPlotOptions({ ...plotOptions, sizeBy: null });
    }
  }, [availableProperties, availableOpacityProperties, plotOptions]);

  React.useEffect(() => {
    if (!collections.find((c) => c.name === "ChEMBL")) {
      // TODO: Add prop to manage additional points
      getChemblUMAPEmbedding()
        .then((data) => setCollections([...collections, { data, name: "ChEMBL" }]))
        .catch((e) => console.error("Error loading chembl umap", e));
    }
  }, [collections, setCollections]);

  const additionalTracesFunction = React.useCallback(
    (particles: IParticle[], get: (p: IParticle | null, axis: "x" | "y") => number) => {
      void clusters; // Line only exists to avoid eslint-disable-next-line react-hooks/exhaustive-deps
      return null;
      // if (!connectByValues || connectByValues.length === 0) {
      //   return [];
      // }

      // const selectedInstancesByConnectBy = connectByValues.map(
      //   (connectBy) =>
      //     new Set([...selected, hover].map((p) => p?.properties?.[connectBy]?.toString()).filter((id) => id != null))
      // );

      // const showHoverOnly: boolean = true;

      // const filteredParticles = showHoverOnly
      //   ? particles.filter((p) =>
      //       connectByValues.every((connectBy, i) =>
      //         selectedInstancesByConnectBy[i].has(p.properties?.[connectBy]?.toString())
      //       )
      //     )
      //   : particles;

      // // TOOD: This could be memoized.
      // const groups = Object.entries(
      //   groupByLodash<IParticle>(filteredParticles, (p) =>
      //     connectByValues.map((connectBy) => `${connectBy}:${p.properties?.[connectBy]}`).join(", ")
      //   )
      // );

      // const allInstances: (IParticle | null)[] = [];
      // const allSizes: (number | null)[] = [];
      // const allColors: (string | null)[] = [];

      // const lineOpacityScaling = scaleSymlog().range([0.1, 0.8]).domain([particles.length, 0]);
      // const hoverColor = color(/* hover?.plotData?.color ||  */ "rgb(0,0,0)")!;

      // const lineColor = hoverColor.copy();
      // lineColor.opacity = lineOpacityScaling(filteredParticles.length);
      // const markerBorderColor = "gray";
      // // hoverColor.opacity = 0.5;
      // const sizeScaling = scaleLinear().range([6, 12]);
      // // Cool plotly optimization: instead of creating many traces for lines, create a single trace with NaN separators.
      // // See https://www.somesolvedproblems.com/2019/05/how-to-make-plotly-faster-with-many.html
      // for (let [key, instances] of groups) {
      //   instances = instances.filter((p) => get(p, "x") != null && get(p, "y") != null);
      //   const instanceScaler = sizeScaling.domain([0, instances.length]);
      //   allInstances.push(...instances);
      //   allInstances.push(null);
      //   allSizes.push(...instances.map((_, i) => instanceScaler(i)));
      //   allSizes.push(null);
      //   allColors.push(...instances.map((p, i) => (p === hover || p.selected ? "gold" : "lightgray")));
      //   allColors.push("darkblue");
      // }

      // return [
      //   {
      //     type: "scattergl",
      //     mode: "lines+markers", // TODO: Maybe use lines+markers
      //     x: allInstances.map((p) => get(p, "x") ?? NaN),
      //     y: allInstances.map((p) => get(p, "y") ?? NaN),
      //     hoverinfo: "all",
      //     marker: {
      //       color: allColors,
      //       size: allSizes,
      //       line: {
      //         color: markerBorderColor.toString(),
      //         width: 2,
      //       },
      //       opacity: 1,
      //     },
      //     line: {
      //       color: lineColor.toString(),
      //       // opacity: 0.5
      //       // shape: "spline",
      //     },
      //     showlegend: false,
      //   },
      // ] as Plotly.Data[];

      // // const result: Plotly.Data[] = [];
      // // groups.forEach(([group, instances]) => {
      // //   // const c = color(
      // //   //   instances.find((i) => i.plotData?.color)?.plotData?.color || "#000000"
      // //   // )!;
      // //   const c = color(hover?.plotData?.color || "#000000")!;
      // //   const sizeScaling = scaleLinear()
      // //     .domain([0, instances.length])
      // //     .range([8, 16]);
      // //   c.opacity = hover && instances.includes(hover) ? 0.5 : 0.3;
      // //   if (true || instances.length > 20 || groups.length > 20) {
      // //     result.push({
      // //       type: "scattergl",
      // //       mode: "lines+markers",
      // //       x: instances.map((p) => get(p, "x")),
      // //       y: instances.map((p) => get(p, "y")),
      // //       name: group,
      // //       hoverinfo: "all",
      // //       marker: {
      // //         color: "lightblue",
      // //         // @ts-ignore
      // //         size: instances.map((p, i) => sizeScaling(i)),
      // //         line: {
      // //           color: c.toString(),
      // //           width: 2,
      // //         },
      // //         opacity: 1,
      // //         // line: 'green'
      // //       },
      // //       line: {
      // //         color: c.toString(),
      // //         // shape: "spline",
      // //       },
      // //       showlegend: false,
      // //     });
      // //   } else {
      // //     // const widthScaling = scaleLinear()
      // //     //   .domain([0, instances.length])
      // //     //   .range([1, 6]);
      // //     // instances.slice(1).forEach((instance, i) => {
      // //     //   const instancePair = [instances[i], instance];
      // //     //   result.push({
      // //     //     type: "scattergl",
      // //     //     mode: "lines+markers",
      // //     //     x: instancePair.map((p) => get(p, "x")),
      // //     //     y: instancePair.map((p) => get(p, "y")),
      // //     //     name: group,
      // //     //     hoverinfo: "all",
      // //     //     marker: {
      // //     //       color: 'red',
      // //     //       size: 10
      // //     //     },
      // //     //     line: {
      // //     //       width: widthScaling(i),
      // //     //       color: c.toString(),
      // //     //       shape: "spline",
      // //     //       // dash: 'dot',
      // //     //     },
      // //     //     // legendgroup: `${group}`,
      // //     //     // showlegend: i === 0,
      // //     //     showlegend: false,
      // //     //   });
      // //     // });
      // //   }
      // // });

      // // return result;
    },
    [
      plotOptions.connectByValues,
      selection,
      clusters,
      // TODO: Hover is excluded to avoid rerendering. Find better solution.
      // hover
    ]
  );

  // const visibleStructures = React.useMemo<IParticle[]>(
  //   () => (showSelectedOnly && selected.length > 0 ? mainParticles.filter((p) => selected.includes(p)) : mainParticles),
  //   [showSelectedOnly, selected, mainParticles]
  // );

  const filteredCollections = React.useMemo<ICollection[]>(() => {
    return filtered ? Object.entries(filtered).map(([name, data]) => ({ name, data })) : collections;
  }, [filtered, collections]);

  const visibleCollections = React.useMemo<ICollection[]>(() => {
    return showSelectedOnly && selection
      ? Object.entries(selection).map(([name, data]) => ({ name, data }))
      : collections;
  }, [collections, showSelectedOnly, selection]);

  return (
    <>
      <Tooltip>
        {hover ? (
          <StructureCard
            structure={hover}
            showProperties={false}
            style={{
              width: "180px",
            }}
          />
        ) : null}
      </Tooltip>
      <HorizontalCollapse
        label="Options"
        position="left"
        size="col-md-2"
        collapsed={optionsCollapsed}
        setCollapsed={setOptionsCollapsed}
      >
        <ComputeEmbeddingsForm
          addCollection={(collection) => setCollections([...collections, collection])}
          loading={loading}
          setLoading={setLoading}
        />
        <MSOForm
          addCollection={(collection) => setCollections([...collections, collection])}
          selection={selection}
          loading={loading}
          setLoading={setLoading}
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
        {/* <RecomputeEmbeddingsForm
          interpolatedStructures={interpolatedStructures}
          setInterpolatedStructures={setInterpolatedStructures}
          particles={mainParticles}
          setParticles={setMainCollection}
          loading={loading}
          setLoading={setLoading}
        /> */}
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
      <div
        // className="col-md-10"
        style={{
          display: "flex",
          flex: 1,
          flexDirection: "column",
          marginLeft: 33,
          marginRight: 33,
          height: "100%",
          overflow: "auto",
          // overflow
        }}
      >
        <LoadingPage
          loading={loading}
          fallback="Please select structures for embedding"
          loadingText="Computing embedding of structures..."
        >
          <details>
            <summary>Projections ({enabledProjections.length})</summary>
            <form className="form-inline" style={{ alignItems: "baseline", flexFlow: "row" }}>
              {availableProjections.map((projection) => (
                <div className="form-check form-check-inline">
                  <input
                    className="form-check-input"
                    type="checkbox"
                    id={`projection${projection}Checkbox`}
                    checked={enabledProjections.includes(projection)}
                    onClick={(e) =>
                      setEnabledProjections(
                        e.currentTarget.checked
                          ? [...enabledProjections, projection]
                          : enabledProjections.filter((p) => p !== projection)
                      )
                    }
                  />
                  <label className="form-check-label" htmlFor={`projection${projection}Checkbox`}>
                    {projection}
                  </label>
                </div>
              ))}
            </form>
          </details>
          {collections.length > 0 ? (
            <>
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
                    {name} Options ({Object.values(plotOptions).filter(Boolean).length})
                  </summary>
                  <form className="form-inline" style={{ alignItems: "baseline", flexFlow: "row" }}>
                    <PlotSelect
                      id={`colorEncoding${i}Select`}
                      label="Color Encoding"
                      // disabled={groupBy !== null}
                      option={plotOptions.colorBy}
                      options={availableProperties}
                      setOption={(colorCoding: string) => {
                        setter({ colorBy: colorCoding });
                      }}
                    />
                    <PlotSelect
                      id={`opacity${i}Select`}
                      label="Opacity By"
                      // disabled={colorCoding !== null}
                      option={typeof plotOptions.opacityBy === "number" ? "constant" : plotOptions.opacityBy}
                      options={availableOpacityProperties}
                      setOption={(opacityBy: string) => {
                        setter({ opacityBy: opacityBy === "constant" ? 0.5 : opacityBy });
                      }}
                    />
                    {typeof plotOptions.opacityBy === "number" ? (
                      <div className="form-group mr-sm-2">
                        {/* <label for="formControlRange">Example Range input</label> */}
                        <input
                          type="range"
                          className="form-control-range"
                          min={0}
                          max={1}
                          step={0.01}
                          value={plotOptions.opacityBy}
                          onChange={(e) => {
                            setter({ opacityBy: e.currentTarget.valueAsNumber });
                          }}
                        />
                      </div>
                    ) : null}
                    <PlotSelect
                      id={`grouping${i}Select`}
                      label="Group By"
                      // disabled={colorCoding !== null}
                      option={plotOptions.groupBy}
                      options={availableProperties}
                      setOption={(groupBy: string) => {
                        setter({ groupBy });
                      }}
                    />
                    <PlotSelect
                      id={`sizeBy${i}Select`}
                      label="Size by"
                      option={plotOptions.sizeBy}
                      options={availableProperties}
                      setOption={(sizeBy: string) => {
                        setter({ sizeBy });
                      }}
                    />
                    <PlotSelect
                      id={`connectBy${i}Select`}
                      label="Connect by"
                      multi={true}
                      option={plotOptions.connectByValues}
                      options={availableProperties}
                      setOption={(option: string[]) => {
                        setter({ connectByValues: option.length === 0 ? null : option });
                      }}
                    />
                  </form>
                </details>
              ))}
              <div className="row m-0" style={{ position: "relative" }}>
                <Grid>
                  <GridItemOptions
                    key="Custom"
                    gridOptions={{
                      w: 12,
                      h: 25,
                      y: 0,
                    }}
                  >
                    <ScatterPlot
                      title={customX && customY ? `${customX} vs. ${customY}` : "Please select two axis properties"}
                      collections={filteredCollections}
                      options={plotOptions}
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
                        <PlotSelect
                          options={availableProperties}
                          option={customY}
                          setOption={(v: string) => setCustomY(v)}
                          id="customYSelect"
                        />
                      </div>
                      <div
                        style={{
                          width: 150,
                          position: "absolute",
                          left: "40%",
                          bottom: 0,
                        }}
                      >
                        <PlotSelect
                          options={availableProperties}
                          option={customX}
                          setOption={(v: string) => setCustomX(v)}
                          id="customXSelect"
                        />
                      </div>
                    </ScatterPlot>
                  </GridItemOptions>
                  {
                    (enabledProjections.length > 0 ? enabledProjections : availableProjections).map((projection, i) => (
                      <GridItemOptions
                        key={projection}
                        gridOptions={{
                          w: 6,
                          h: 25,
                          y: 25 + Math.floor(i / 2) * 25,
                          x: i % 2 === 0 ? 0 : 6,
                        }}
                      >
                        <ScatterPlot
                          title={projection}
                          collections={filteredCollections}
                          options={plotOptions}
                          setHover={setHover}
                          selected={selection}
                          setSelected={setSelection}
                          xAccessor={`projection[${projection}][0]`}
                          yAccessor={`projection[${projection}][1]`}
                        />
                      </GridItemOptions>
                    )) as any
                  }
                  {/* <GridItemOptions
                    key="TMAP"
                    gridOptions={{
                      w: 6,
                      h: 25,
                      y: 25,
                      x: 0,
                    }}
                  >
                    <ProjectionPlot
                      title="TMAP"
                      markerFunction={markerFunction}
                      transformsFunction={transformsFunction}
                      particles={filtered || mainParticles}
                      setHover={setHover}
                      selected={selected}
                      setSelected={setSelected}
                      additionalParticles={interpolatedStructures}
                      additionalTracesFunction={additionalTracesFunction}
                      xAccessor="projection[tmap][0]"
                      yAccessor="projection[tmap][1]"
                    />
                  </GridItemOptions>
                  <GridItemOptions
                    key="UMAP"
                    gridOptions={{
                      w: 6,
                      h: 25,
                      y: 25,
                      x: 6,
                    }}
                  >
                    <ProjectionPlot
                      title="UMAP"
                      markerFunction={markerFunction}
                      transformsFunction={transformsFunction}
                      particles={filtered || mainParticles}
                      setHover={setHover}
                      selected={selected}
                      setSelected={setSelected}
                      additionalParticles={interpolatedStructures}
                      additionalTracesFunction={additionalTracesFunction}
                      xAccessor="projection[umap][0]"
                      yAccessor="projection[umap][1]"
                    />
                  </GridItemOptions>
                  <GridItemOptions
                    key="ChEMBL UMAP"
                    gridOptions={{
                      w: 6,
                      h: 25,
                      y: 25,
                      x: 6,
                    }}
                  >
                    <UMAPProjectionPlot
                      umapPoints={umapPoints}
                      colorCoding={colorCoding}
                      markerFunction={markerFunction}
                      transformsFunction={transformsFunction}
                      particles={filtered || mainParticles}
                      setHover={setHover}
                      selected={selected}
                      setSelected={setSelected}
                      additionalParticles={interpolatedStructures}
                      additionalTracesFunction={additionalTracesFunction}
                      xAccessor="projection[chembl_umap][0]"
                      yAccessor="projection[chembl_umap][1]"
                    />
                  </GridItemOptions> */}
                  {/* <GridItemOptions
                    key="ChEMBL TSNE"
                    gridOptions={{
                      w: 6,
                      h: 25,
                      y: 25,
                      x: 6,
                    }}
                  >
                    <ProjectionPlot
                      title="Fixed TSNE"
                      markerFunction={markerFunction}
                      transformsFunction={transformsFunction}
                      particles={filtered || particles}
                      setHover={setHover}
                      selected={selected}
                      setSelected={setSelected}
                      additionalParticles={interpolatedStructures}
                      additionalTracesFunction={additionalTracesFunction}
                      xAccessor="projection[chembl_tsne][0]"
                      yAccessor="projection[chembl_tsne][1]"
                    />
                  </GridItemOptions>
                  <GridItemOptions
                    key="TSNE"
                    gridOptions={{
                      w: 6,
                      h: 25,
                      y: 50,
                      x: 0,
                    }}
                  >
                    <ProjectionPlot
                      title="TSNE"
                      markerFunction={markerFunction}
                      transformsFunction={transformsFunction}
                      particles={filtered || particles}
                      setHover={setHover}
                      selected={selected}
                      setSelected={setSelected}
                      additionalParticles={interpolatedStructures}
                      additionalTracesFunction={additionalTracesFunction}
                      xAccessor="projection[tsne][0]"
                      yAccessor="projection[tsne][1]"
                    />
                  </GridItemOptions>
                  <GridItemOptions
                    key="PCA"
                    gridOptions={{
                      w: 6,
                      h: 25,
                      y: 50,
                      x: 6,
                    }}
                  >
                    <ProjectionPlot
                      title="PCA"
                      markerFunction={markerFunction}
                      transformsFunction={transformsFunction}
                      particles={filtered || particles}
                      setHover={setHover}
                      selected={selected}
                      setSelected={setSelected}
                      additionalParticles={interpolatedStructures}
                      additionalTracesFunction={additionalTracesFunction}
                      xAccessor="projection[pca][0]"
                      yAccessor="projection[pca][1]"
                    />
                  </GridItemOptions> */}
                </Grid>
              </div>
              <StructureCardGrid
                collections={visibleCollections}
                selected={selection}
                setSelected={setSelection}
                setFiltered={setFiltered}
                clusters={clusters}
                tableClass="main-ranking"
                structureCardProps={(structure) => ({
                  className: structure === hover ? "border-primary" : structure.selected ? "border-secondary" : "",
                })}
                renderTopForm={
                  <>
                    <div className="custom-control custom-switch mr-sm-2">
                      <input
                        type="checkbox"
                        className="custom-control-input"
                        id="customSwitch1"
                        checked={showSelectedOnly}
                        onChange={(e) => setShowSelectedOnly(e.currentTarget.checked)}
                      />
                      <label className="custom-control-label" htmlFor="customSwitch1">
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
      <HorizontalCollapse
        label="Clusters"
        position="right"
        size="col-md-2"
        collapsed={clusterCollapsed}
        setCollapsed={setClusterCollapsed}
      >
        <ClusterSidePanel
          clusters={clusters}
          setClusters={setClusters}
          selected={selection}
          setSelected={setSelection}
        />
      </HorizontalCollapse>
    </>
  );
}
