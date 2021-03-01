import * as React from "react";
import LineUp, {
  builder,
  buildRanking,
  buildStringColumn,
  Column,
  equal,
  IColumnDesc,
  LocalDataProvider,
  Taggle,
  buildCategoricalColumn,
  Ranking,
  StringColumn,
  IStringFilter,
} from "lineupjs";
import { ICollection, IParticle, IParticleSelection } from "../../interfaces";
import debounce from "lodash.debounce";
import { useSyncedRef } from "../../utils/hooks";
import { StructureImageRenderer } from "./StructureImageRenderer";
import { getClustersFromParticle } from "../../ClusterSidePanel";
import { StructureImageColumn } from "./StructureImageColumn";
import { hasSubstructureMatch } from "../../utils/api";

export const LineupWrapper = React.memo(
  ({
    collections,
    setSelected,
    selected,
    setFiltered,
    clusters,
    ...innerProps
  }: {
    collections: ICollection[];
    setSelected?: (selection: IParticleSelection) => void;
    selected?: IParticleSelection;
    setFiltered?: (filtered: IParticleSelection) => void;
    clusters?: { [key: string]: IParticle[] };
  } & React.HTMLAttributes<HTMLDivElement>) => {
    const divRef = React.useRef<HTMLDivElement>(null);
    const lineupRef = React.useRef<Taggle | null>(null);
    const rankingRef = React.useRef<Ranking | null>(null);

    const setSelectedRef = useSyncedRef(setSelected);
    const setFilteredRef = useSyncedRef(setFiltered);

    // Reduce the collections to an array of objects with a dataset property
    const mergedData = React.useMemo(
      () =>
        collections.reduce<(IParticle & { _dataset: string; _particle: IParticle })[]>((acc, cur) => {
          cur.data.forEach((d) =>
            acc.push({
              ...d,
              _dataset: cur.name,
              _particle: d,
            })
          );
          return acc;
        }, []),
      [collections]
    );

    React.useEffect(() => {
      /* const previousDataSet = new Set(previousData);
      if(lineupRef.current && data.some((d) => previousDataSet.has(d))) {
        console.log("SHALLOW");
        rankingRef.current?.markDirty('all');
        lineupRef.current?.update();
      } else  if(!lineupRef.current || previousData !== data) {
        */
      lineupRef.current?.destroy();

      const DEFAULT_HEIGHT = 18;
      let height = 36;

      const b = builder(mergedData)
        .ranking(buildRanking().supportTypes().allColumns())
        .aggregationStrategy("group+item+top")
        .propagateAggregationState(true)
        .registerColumnType("structureImage", StructureImageColumn)
        .registerRenderer("structureImage", new StructureImageRenderer())
        .livePreviews({
          filter: false,
        })
        .deriveColumns("_dataset")
        .column(
          buildCategoricalColumn("properties.clusters")
            .asSet(";")
            .categories(Array.from(new Set(mergedData.map((p) => getClustersFromParticle(p)).flat())))
            .label("Clusters")
            .build([])
        )
        // .deriveColumns("structure")
        .column(
          buildStringColumn("structure")
            .label("Structure")
            // eslint-disable-next-line no-template-curly-in-string
            .pattern("https://pubchem.ncbi.nlm.nih.gov/#query=${value}")
        )
        .column({
          ...buildStringColumn("structure")
            .label("Structure")
            .renderer("structureImage", "structureImage")
            .width(height * 3)
            .build([]),
          type: "structureImage",
        })
        .sidePanel(true, true)
        .deriveColors()
        .dynamicHeight(() => ({
          defaultHeight: DEFAULT_HEIGHT,
          height: () => height,
          padding: () => 0,
        }));

      // TODO: I really don't like that, but Lineup infers the columns by the *first* entry in the data.
      // So here we move the object with the most properties to the first position. We should probably create a "merged" object for that.
      // https://github.com/lineupjs/lineupjs/blob/develop/src/provider/utils.ts#L268-L269
      const mergedDataProperties = mergedData.map(({ properties = {} }) => properties);
      const propertyWithMostEntries =
        mergedDataProperties.reduce<object | null>((acc, cur) => {
          return !acc || Object.keys(acc!).length < Object.keys(cur).length ? cur : acc;
        }, null) || {};

      // The properties are a nested object, so we derive it using the builder and then inject the 'properties.' in the column.
      const propertiesBuilder = builder([propertyWithMostEntries, ...mergedDataProperties]).deriveColumns();

      // @ts-ignore
      const propertiesColumns = (propertiesBuilder.columns as IColumnDesc[]).map(
        // @ts-ignore
        (col) => ({ ...col, column: `properties.${col.column}` })
      );
      propertiesColumns.filter((col) => col.label !== "Clusters").forEach((col) => b.column(col));

      // Change the default renderer of the embedding column
      // b.deriveColumns("embedding");
      // @ts-ignore
      (b.columns as IColumnDesc[]).forEach((col) => {
        // Patch the default renderer for all numbers columns
        if (col.type === "numbers") {
          col.renderer = "histogram"; // verticalbar is nice but slow..
          col.groupRenderer = "histogram"; // verticalbar is nice but slow..
        }
      });

      // Build the ranking
      const lineup = b.buildTaggle(divRef.current!);

      // Lookup the image column
      const imageColumn: StructureImageColumn | null = lineup.data.find(
        (col) => col.getRenderer() === "structureImage"
      ) as StructureImageColumn | null;
      if (!imageColumn) {
        console.error("Column with image renderer not found. Autoresizing disabled.");
      }
      const debouncedWidthChanged = debounce((prev, cur) => {
        // Adjust the height according to the image column
        height = Math.min(Math.max(DEFAULT_HEIGHT, (cur - 0) / 3), 150);
        lineup.update();
      }, 500);

      const debouncedFilterChanged = debounce((prev, cur: IStringFilter) => {
        const filter = typeof cur?.filter === "string" ? cur?.filter : null;
        if (imageColumn && filter) {
          hasSubstructureMatch(
            mergedData.map((d) => d.structure),
            filter
          )
            .then((matches) => {
              // TODO: Race condition check
              const validSmiles = Object.entries(matches.counts)
                .filter(([smiles, count]) => count > 0)
                .map(([smiles, count]) => smiles);
              imageColumn.setFilter({
                filter,
                valid: new Set(validSmiles),
                filterMissing: true,
              });
            })
            .catch((e) => {
              console.log(e);
              imageColumn.setFilter(null);
            });
        }
      }, 1000);
      imageColumn?.on(Column.EVENT_WIDTH_CHANGED, debouncedWidthChanged);
      imageColumn?.on(StringColumn.EVENT_FILTER_CHANGED, debouncedFilterChanged);

      const rowsToSelection = (rows: any[]) => {
        return rows.reduce<{ [key: string]: IParticle[] }>((acc, cur) => {
          if (!acc[cur._dataset]) {
            acc[cur._dataset] = [];
          }
          acc[cur._dataset].push(cur._particle);
          return acc;
        }, {});
      };

      // Listen to selection
      lineup.on(LineUp.EVENT_SELECTION_CHANGED, async () => {
        const data = rowsToSelection(await lineup.data.view(lineup.getSelection()));
        setSelectedRef.current?.(Object.entries(data).length === 0 ? null : data);
      });

      lineup.data.on(LocalDataProvider.EVENT_ORDER_CHANGED, async (oldSelection, newSelection) => {
        if (newSelection.length === mergedData.length) {
          setFilteredRef?.current?.(null);
        } else if (!equal(oldSelection.sort(), newSelection.sort())) {
          // If we actually filtered the table
          setFilteredRef?.current?.(rowsToSelection(await lineup.data.view(newSelection)));
        }
      });

      // @ts-ignore
      lineupRef.current = lineup;
      rankingRef.current = lineup.data.getFirstRanking();
      // }

      return () => {
        // lineupRef.current?.destroy();
      };
    }, [collections, mergedData, clusters, setFilteredRef, setSelectedRef]);

    return (
      <div
        {...(innerProps || {})}
        ref={divRef}
        style={{
          ...(innerProps?.style || {}),
          flex: 1,
          fontSize: "smaller",
        }}
      ></div>
    );
  }
);
