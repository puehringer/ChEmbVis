import * as React from "react";
import { Figure } from "react-plotly.js";
import { ICollection, IParticleSelection } from "../interfaces";
import { PlotComponent, PLOTLY_CONFIG } from "./PlotComponent";
import { max, min } from "d3-array";
import Select, { OptionsType } from "react-select";
import isEqual from "lodash.isequal";
import { GridItemOptions } from "./GridItemOptions";
import { Grid } from "./Grid";
import { DEFAULT_COLORWAY, DEFAULT_PLOTLY_COLORSCALE } from "../utils/constants";

export interface IParallelCoordinatesPlotProps {
  collections: ICollection[];
  selection: IParticleSelection;
  setSelection(selected: IParticleSelection): void;
  /**
   * Currently filtered selection. If passed, only the coordinates of filtered items are shown.
   */
  filterSelection?: IParticleSelection;
}

interface ISelectOption {
  label: string;
  value: string;
  order: "asc" | "desc";
}

export const ParallelCoordinatesPlot = React.memo(({ collections, selection, setSelection, filterSelection }: IParallelCoordinatesPlotProps) => {
  const [figureState, setFigureState] = React.useState<Figure | null>(null);
  const [enabledProperties, setEnabledProperties] = React.useState<OptionsType<ISelectOption>>([]);
  const [colorProperty, setColorProperty] = React.useState<ISelectOption | null>(null);
  // const [selection, setSelection] = React.useState<IParticleSelection | null>(null);
  const [constraintRange, setConstraintRange] = React.useState<{ [key: string]: [number, number][] } | null>(null);

  const filteredCollections = React.useMemo(() => collections.filter((c) => !c.hidden), [collections]);
  const availableProperties = React.useMemo(
    () =>
      Array.from(
        new Set(
          filteredCollections
            .map((c) =>
              Object.keys(c.data[0]?.properties || {})
                // TODO: Find other condition
                // .filter(([key, value]) => typeof value === "number")
                // .map(([key, value]) => key)
            )
            .flat()
        )
      ),
    [filteredCollections]
  );

  React.useEffect(() => {
    // Cleanup
    // setFilterSelection(null);
    setConstraintRange(null);
    // setSelection(null);
  }, [filteredCollections, enabledProperties]);

  React.useEffect(() => {
    setFigureState((figureState) => {
      if (enabledProperties.length < 2) {
        return null;
      }

      const filtered = Object.values(filterSelection || {}).flat();
      const all = filtered.length > 0 ? filtered : filteredCollections.map((c) => c.data).flat().slice(0, 500);

      const newFigureState: Figure = {
        frames: [],
        layout: {
          ...(figureState?.layout || {}),
          autosize: true,
          colorway: DEFAULT_COLORWAY,
          showlegend: true,
        },
        data: [
          {
            // @ts-ignore
            type: "parcoords",
            line: {
              //   // TODO Make semi-transparent,
              // @ts-ignore
              showscale: colorProperty ? false : false,
              // colorscale: [
              //   ['0.0', DEFAULT_COLORWAY[0] + '0a'],
              //   ['1.0', DEFAULT_COLORWAY[1] + '0a']
              // ],
              // If we don't have a numeric scale, use a color map for the filteredCollections instead
              // colorscale: colorProperty ? "Portland" : filteredCollections.map((c, i) => [i, DEFAULT_COLORWAY[i]]),
              ...DEFAULT_PLOTLY_COLORSCALE,
              color: colorProperty ? all.map((p) => p.properties![colorProperty.value] as number) : filteredCollections.map((c, i) => Array(c.data.length).fill(i.toString())).flat(),
            },
            // @ts-ignore
            dimensions: [
              ...enabledProperties.map(({ value, order }) => {
                const values = all.map((p) => p.properties![value] as number);
                const minValue = min(values);
                const maxValue = max(values);
                if (minValue == null || maxValue == null) {
                  return null;
                }
                const range = [Math.min(0, minValue), maxValue];

                return {
                  range: order === "desc" ? range : range.reverse(),
                  constraintrange: constraintRange?.[value],
                  values,
                  label: value,
                };
              }),
            ],
          },
        ],
      };
      console.log(newFigureState)
      return newFigureState;
    });
  }, [filteredCollections, filterSelection, enabledProperties, constraintRange, colorProperty]);

  return (
    <div className="d-flex flex-column flex-fill">
      <form>
        <div className="row">
          <div className="col-md-10 mb-3">
            <label htmlFor="propertiesSelect">Properties</label>
            <Select<ISelectOption, true>
              menuPosition="fixed"
              isMulti
              className="mb-3"
              name="propertiesSelect"
              value={enabledProperties}
              onChange={(e) => {
                setEnabledProperties(e);
              }}
              options={availableProperties.map((p) => ({
                value: p,
                label: p,
                order: "desc",
              }))}
              openMenuOnClick={false}
              closeMenuOnSelect={false}
              formatOptionLabel={(option, meta) => {
                if (meta.context === "value") {
                  return (
                    <div className="d-flex flex-row align-items-center">
                      <span title={option.label} className="text-truncate" style={{ maxWidth: 150 }}>
                        {option.label}
                      </span>
                      <i
                        role="button"
                        title="Change axis order"
                        className={`fas fa-fw fa-sort-amount-${option.order === "asc" ? "up-alt" : "down"} ms-1 me-1 `}
                        onClick={(e) => {
                          setEnabledProperties(
                            enabledProperties.map((p) =>
                              p === option
                                ? {
                                    ...p,
                                    order: p.order === "asc" ? "desc" : "asc",
                                  }
                                : p
                            )
                          );
                        }}
                      />
                    </div>
                  );
                }
                return option.label;
              }}
            />
          </div>
          <div className="col-md-2 mb-3">
            <label htmlFor="colorBySelect">Color By</label>
            <Select<ISelectOption, false>
              menuPosition="fixed"
              className="mb-3"
              name="colorBySelect"
              value={colorProperty}
              isClearable={true}
              onChange={(e) => {
                setColorProperty(e);
              }}
              options={availableProperties.map((p) => ({
                value: p,
                label: p,
                order: "desc",
              }))}
            />
          </div>
        </div>
      </form>

      {figureState ? (
        <>
          <Grid>
            <GridItemOptions
              enableMove={false}
              key="parcoordsGridItem"
              gridOptions={{
                w: 12,
                h: 25,
                y: 0,
                x: 0,
              }}
            >
            <i
              className="fas fa-fw fa-sync-alt react-grid-item-hidden"
              title="Sync global selection with plot (up to 1000)"
              onClick={() => {
                const allSelected = Object.values(selection || {}).flat();
                if (allSelected.length <= 1000) {
                  // If we have below a certain threshold of items, we can create the constraint range for each individual one
                  setConstraintRange(
                    enabledProperties.reduce((acc, { value }) => {
                      return {
                        ...acc,
                        [value]: allSelected.map((p) => p.properties![value] as number).map((v) => [v, v + 0.0001]),
                      };
                    }, {})
                  );
                } else {
                  setConstraintRange(null);
                  // TODO: Check if current range is sufficient for current selection
                  // setConstraintRange((currentRange) => {
                  //   const allIncluded = allSelected.every((p) => Object.entries(currentRange || {}).every(([key, allRanges]) => {
                  //     return allRanges.some(
                  //       (range) => p.properties?.[key]! >= range[0] && p.properties?.[key]! <= range[1]
                  //     );
                  //   }));
                  //   return allIncluded ? currentRange : null;
                  // });
                }
              }}
              style={{
                position: "absolute",
                zIndex: 1,
                top: 3,
                left: 55,
                cursor: "pointer",
              }}
            />
            <i
              className="fas fa-fw fa-ban react-grid-item-hidden"
              title="Clear selection"
              onClick={() => {
                setConstraintRange(null);
              }}
              style={{
                position: "absolute",
                zIndex: 1,
                top: 3,
                left: 80,
                cursor: "pointer",
              }}
            />
              <PlotComponent
                style={{
                  width: "99%",
                  height: "100%",
                }}
                data={figureState.data}
                layout={figureState.layout}
                config={PLOTLY_CONFIG}
                onRestyle={() => {
                  // @ts-ignore
                  const filter: { [key: string]: [number, number][] } = figureState.data[0].dimensions
                    // @ts-ignore
                    .filter((dim) => dim?.constraintrange)
                    .reduce(
                      // @ts-ignore
                      (acc, cur) =>
                        cur.constraintrange?.length > 0 ? { ...acc, [cur.label]: typeof cur.constraintrange[0] === 'number' ? [cur.constraintrange] : cur.constraintrange } : acc,
                      {}
                    );

                  if (Object.entries(filter).length > 0) {
                    const valid = filteredCollections.reduce(
                      (acc, c) => ({
                        ...acc,
                        [c.name]: c.data.filter((p) =>
                          Object.entries(filter).every(([key, allRanges]) => {
                            return allRanges.some(
                              (range) => p.properties?.[key]! >= range[0] && p.properties?.[key]! <= range[1]
                            );
                          })
                        ),
                      }),
                      {}
                    );
                    setConstraintRange(filter);
                    setSelection(Object.entries(valid).length > 0 ? valid : null);
                  } else {
                    setConstraintRange(null);
                    setSelection(null);
                  }
                }}
                onUpdate={(figure) => {
                  // Synchronize ordering of parallel coordinates with select
                  // @ts-ignore
                  const newOrdering: string[] = figure.data?.[0].dimensions.filter(Boolean).map((d) => d.label);
                  if (
                    newOrdering &&
                    !isEqual(
                      newOrdering,
                      enabledProperties.map((e) => e.value)
                    )
                  ) {
                    setEnabledProperties(
                      newOrdering.map((value) => enabledProperties.find((e) => e.value === value)!).filter(Boolean)
                    );
                  }
                  // Inline save the figure to save the zoom
                  //   figureState.data = figure.data;
                  //   figureState.frames = figure.frames;
                  //   figureState.layout = figure.layout;
                }}
              />
            </GridItemOptions>
          </Grid>
        </>
      ) : (
        <p>Please select at least 2 properties from the dropdown.</p>
      )}
    </div>
  );
});
