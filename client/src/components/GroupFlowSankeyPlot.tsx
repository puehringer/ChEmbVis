import * as React from "react";
import Select from "react-select/creatable";
import { ICollection, IParticle, IParticleSelection } from "../interfaces";
import { Tooltip } from "../utils/hooks";
import { PlotComponent, PLOTLY_CONFIG } from "./PlotComponent";
import { LineupWrapper } from "./ranking/LineupWrapper";
import { StructureImage } from "./StructureImage";

const LAYOUT: Partial<Plotly.Layout> = {
  dragmode: "lasso",
  hovermode: "closest",
  autosize: true,
  margin: {
    // l: 250,
    r: 0,
    b: 0,
    t: 50,
    // pad: 4,
  },
};

export const GroupFlowSankeyPlot = React.memo(({
  collections,
  selection,
}: //   setHover,
//   hover,
{
  selection: IParticleSelection;
  collections: ICollection[];
}) => {
  const [collection, setCollection] = React.useState<ICollection | null>(null);
  const [enabledProperties, setEnabledProperties] = React.useState<string[]>([]);
  const [selectionCollections, setSelectionCollections] = React.useState<ICollection[]>([]);
  const allAvailableProperties = React.useMemo(() => Object.keys(collection?.data?.[0]?.properties || {}), [collection]);
  const [hover, setHover] = React.useState<IParticle[] | null>(null);

  const data = React.useMemo<Partial<Plotly.PlotData>[][]>(() => {
    
    if (!collection || enabledProperties.length < 2) {
      return [];
    }

    const map = new Map<IParticle, number[]>();
    // Each unique value has to have a unique id, i.e. running index in this case. Otherwise, one could not distinguish the same labels of different properties.
    let valueIndex = 0;
    const labels: string[] = [];
    enabledProperties.forEach((property, i) => {
      const valueLookup = new Map<string, number>();
      collection.data.forEach((d) => {
        const value = `${d.properties[property]}`;
        if (!valueLookup.has(value)) {
          valueLookup.set(value, valueIndex++);
          labels.push(`${value}`);
        }

        if (!map.has(d)) {
          map.set(d, []);
        }
        map.get(d)![i] = valueLookup.get(value)!;
      });

      valueLookup.clear();
    });

    // Add selection (if any)
    const selected = selection?.[collection.name];
    if(selected) {
      selected.forEach((particle) => {
        map.get(particle)?.push(valueIndex);
      });

      // Add to labels, increase value index
      labels.push("Selection");
      valueIndex++;
    }

    const matrix: IParticle[][][] = new Array(valueIndex).fill(null).map(() => []);
    map.forEach((assignments, particle) => {
      assignments.slice(1).forEach((next, currentIndex) => {
        const current = assignments[currentIndex];
        matrix[current][next] = matrix[current][next] || [];
        matrix[current][next].push(particle);
      });
    });

    const link = {
      source: [] as number[],
      target: [] as number[],
      value: [] as number[],
      customdata: [] as IParticle[][]
    };
    matrix.forEach((toValues, from) => {
      toValues.forEach((particles, to) => {
        link.source.push(from);
        link.target.push(to);
        link.value.push(particles.length);
        link.customdata.push(particles);
      });
    });

    console.log(matrix);

    return [[
      {
        type: "sankey",
        orientation: "h",
        arrangement: "fixed",
        // @ts-ignore
        node: {
          // pad: 15,
          // thickness: 30,
          // line: {
          // color: "black",
          // width: 0.5,
          // },
          label: labels,
          // groups: [[1,2,3, 4,5,6,7,8,9]]
          // color: ["blue", "blue", "blue", "blue", "blue", "blue"],
          customdata: matrix.map((to, i) => (to[i] || []).flat())
        },
        link: {
          ...link,
          // color: 'rgba(200, 200, 200, 10)'
        },
      },
    ]];
  }, [selection, collection, enabledProperties]);

  return (
    <div
      className="d-flex flex-column"
      style={{
        position: "relative",
        // overflowY: 'scroll',
        // minHeight: 500
      }}
    >
      <Tooltip>
        {hover ? <StructureImage structure={hover.map((h) => h.structure)}
            style={{
              width: 180,
              backgroundColor: "rgba(255, 255, 255, 0.6)",
            }} /> : null}
      </Tooltip>
      <p className="text-muted"></p>
      <div className="row">
        <div className="col-md-3 mb-3">
          <label>Collection</label>
          <Select<{ label: string; value: ICollection }>
            menuPosition="fixed"
            value={collection ? { label: collection.name, value: collection } : null}
            onChange={(e) => {
              setCollection(e?.value || null);
            }}
            options={collections.map((c) => ({
              label: c.name,
              value: c,
            }))}
          />
        </div>
        <div className="col-md-9 mb-3">
          <label>Group by</label>
          <Select<{ label: string; value: string }, true>
            menuPosition="fixed"
            isDisabled={!collection}
            isMulti
            value={enabledProperties.map((p) => ({
              label: p,
              value: p,
            }))}
            onChange={(e) => {
              setEnabledProperties(e.map((p) => p.value));
            }}
            options={allAvailableProperties.map((property) => ({
              label: property,
              value: property,
            }))}
            closeMenuOnSelect={false}
          />
        </div>
      </div>
      {data.length > 0 ? (
        <> 
          {data.map((d, i) => <div><PlotComponent
            key={i}
            className="mt-3 mb-5"
            style={{ width: "100%", height: 500 }}
            data={d}
            layout={LAYOUT}
            config={PLOTLY_CONFIG}
            onClick={(e) => {
              setSelectionCollections(e.points.map((p) => ({data: p.customdata as unknown as IParticle[], name: 'Selection' })));
            }}
            onHover={(e) => {
              const hoveredPoints = (e.points?.map((p) => (p.customdata as unknown as  IParticle[]) || []).flat() || []);
              setHover(hoveredPoints.length > 0 ? hoveredPoints : null);
            }}
            onUnhover={() => {
              setHover(null);
            }}
          /></div>)}
          <div style={{ flex: "1 1 500px", display: "flex", overflow: "auto" }}>
            {selectionCollections ? (
                  <LineupWrapper
                    collections={selectionCollections}
                    // setSelected={setSelectedFromLineup}
                    // getRankingBuilders={enabledEmbeddings.length === 0 ? undefined : getRankingBuilders}
                    // adjustRankings={enabledEmbeddings.length === 0 ? undefined : adjustRankings}
                  />
                ) : null}
          </div>
        </>
      ) : null}
    </div>
  );
});
