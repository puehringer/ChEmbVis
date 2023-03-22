import * as React from "react";
import LineUp, {
  builder,
  buildRanking,
  buildStringColumn,
  Column,
  IColumnDesc,
  LocalDataProvider,
  Taggle,
  Ranking,
  StringColumn,
  IStringFilter,
  RankingBuilder,
} from "lineupjs";
import { ICollection, IParticle, IParticleSelection } from "../../interfaces";
import debounce from "lodash.debounce";
import { useSyncedRef } from "../../utils/hooks";
import { StructureImageRenderer } from "./StructureImageRenderer";
import { StructureImageColumn } from "./StructureImageColumn";
import { hasSubstructureMatch } from "../../utils/api";
import isEqual from "lodash.isequal";
import castArray from "lodash.castarray";

export const buildDefaultRanking = () => {
  const rankingBuilder = buildRanking();
  rankingBuilder.supportTypes();
  rankingBuilder.allColumns();
  return rankingBuilder;
};

export const LineupWrapper = ({
  collections,
  getRankingBuilders = buildDefaultRanking,
  adjustRankings,
  setSelected,
  align,
  selected,
  setFiltered,
  onSelectionSet,
  ...innerProps
}: {
  collections: ICollection[];
  getRankingBuilders?: () => RankingBuilder | RankingBuilder[];
  align?: string;
  adjustRankings?: (lineup: Taggle, rankings: Ranking[]) => void;
  setSelected?: (selection: IParticleSelection) => void;
  selected?: IParticleSelection;
  setFiltered?: (filtered: IParticleSelection) => void;
  onSelectionSet?: (lineup: Taggle, rankings: Ranking[]) => void;
} & React.HTMLAttributes<HTMLDivElement>) => {
  const divRef = React.useRef<HTMLDivElement>(null);
  const lineupRef = React.useRef<Taggle | null>(null);
  const rankingRef = React.useRef<Ranking[] | null>(null);
  const indexMapRef = React.useRef<Map<IParticle, number> | null>(null);
  const disableLineUpSelectionListener = React.useRef<boolean>(false);

  const setSelectedRef = useSyncedRef(setSelected);
  const setFilteredRef = useSyncedRef(setFiltered);
  const onSelectionSetRef = useSyncedRef(onSelectionSet);

  // Reduce the collections to an array of objects with a dataset property
  const mergedData = React.useMemo(
    () =>
      collections.reduce<(IParticle & { _dataset: string; _particle: IParticle } & Record<string, unknown>)[]>((acc, cur) => {
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
    lineupRef.current?.destroy();

    const DEFAULT_HEIGHT = 18;
    let height = 36;

    const rankingBuilders = getRankingBuilders();

    const b = builder(mergedData).animated(false);
    castArray(rankingBuilders).forEach((builder) => b.ranking(builder));
    b.aggregationStrategy("group+item+top")
      .propagateAggregationState(true)
      .registerColumnType("structureImage", StructureImageColumn)
      .registerRenderer("structureImage", new StructureImageRenderer())
      .livePreviews({
        filter: false,
      })
      .sidePanel(true, true)
      // .deriveColors()
      .dynamicHeight(() => ({
        defaultHeight: DEFAULT_HEIGHT,
        height: () => height,
        padding: () => 0,
      }));

    b.deriveColumns("_dataset");

    // Use image column as link and text instead
    // b.column(
    //   buildStringColumn("structure")
    //     .label("Structure")
    //     // eslint-disable-next-line no-template-curly-in-string
    //     .pattern("https://pubchem.ncbi.nlm.nih.gov/#query=${value}")
    //     .build([])
    // );

    b.column({
      ...buildStringColumn("structure")
        .label("Structure")
        .renderer("structureImage", "structureImage")
        .width(height * 3)
        .build([]),
      type: "structureImage",
    });

    // TODO: I really don't like that, but Lineup infers the columns by the *first* entry in the data.
    // So here we move the object with the most properties to the first position. We should probably create a "merged" object for that.
    // https://github.com/lineupjs/lineupjs/blob/develop/src/provider/utils.ts#L268-L269
    const mergedDataProperties = mergedData.map(({ properties = {} }) => properties);
    const propertyWithMostEntries =
      mergedDataProperties.reduce<Record<string, unknown> | null>((acc, cur) => {
        return !acc || Object.keys(acc!).length < Object.keys(cur).length ? cur : acc;
      }, null) || {};

    const lazyColumns = Object.keys(propertyWithMostEntries).filter(
      (key) => Object.getOwnPropertyDescriptor(propertyWithMostEntries, key)!["get"]
    );
    const eagerColumns = Object.keys(propertyWithMostEntries).filter((key) => !lazyColumns.includes(key));
    // console.log(propertyWithMostEntries, lazyColumns, eagerColumns);
    // The properties are a nested object, so we derive it using the builder and then inject the 'properties.' in the column.
    const propertiesBuilder = builder([propertyWithMostEntries, ...mergedDataProperties]).deriveColumns();
    // lazyColumns.forEach((col) =>
    //   propertiesBuilder.column(buildNumberColumn(col).custom("visible", false).custom("lazy", true).build([]))
    // );

    // @ts-ignore
    const propertiesColumns = (propertiesBuilder.columns as IColumnDesc[]).map(
      // @ts-ignore
      (col) => ({ ...col, column: `properties.${col.column}` })
    );
    propertiesColumns.forEach((col) => b.column(col));

    // Change the default renderer of the embedding column
    // TODO: Adjust to multiple embeddings
    // b.deriveColumns("embedding");
    // @ts-ignore
    (b.columns as IColumnDesc[]).forEach((col) => {
      // Patch the default renderer for all numbers columns
      if (col.type === "numbers") {
        col.renderer = "histogram"; // verticalbar is nice but slow..
        col.groupRenderer = "histogram"; // verticalbar is nice but slow..
      }
      if (col.type === "number") {
        col.groupRenderer = "histogram"; // boxplot is nice but slow..
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
            console.error(e);
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
      if (!disableLineUpSelectionListener.current) {
        const data = rowsToSelection(await lineup.data.view(lineup.getSelection()));
        setSelectedRef.current?.(Object.entries(data).length === 0 ? null : data);
      }
    });

    let disableTrigger = false;
    lineup.data.on(LocalDataProvider.EVENT_ADD_COLUMN, (col, index) => {
      if (disableTrigger) {
        return;
      }
      const columnName = col.desc.label.toLowerCase();
      // @ts-ignore
      if (col.desc.lazy) {
        propertiesBuilder.deriveColumns(columnName);
        // @ts-ignore
        const createdColumnDesc = (propertiesBuilder.columns as IColumnDesc[])
          // @ts-ignore
          .filter((col) => col.column === columnName && !col.lazy)
          .map(
            // @ts-ignore
            (col) => ({ ...col, column: `properties.${col.column}` })
          )?.[0];
        disableTrigger = true;
        const createdColumn = lineup.data.getFirstRanking().insertAfter(lineup.data.create(createdColumnDesc)!, col);
        disableTrigger = false;
        createdColumn?.markDirty("all");
        lineup.update();
      }
    });

    lineup.data.on(LocalDataProvider.EVENT_ORDER_CHANGED, async (oldSelection, newSelection) => {
      if (newSelection.length === mergedData.length) {
        setFilteredRef?.current?.(null);
      } else if (!isEqual(oldSelection.slice().sort(), newSelection.slice().sort())) {
        // If we actually filtered the table
        setFilteredRef?.current?.(rowsToSelection(await lineup.data.view(newSelection)));
      }
    });

    // @ts-ignore
    lineupRef.current = lineup;
    rankingRef.current = lineup.data.getRankings();

    adjustRankings?.(lineupRef.current, rankingRef.current);

    indexMapRef.current = mergedData.reduce((acc, cur, i) => {
      acc.set(cur._particle, i);
      return acc;
    }, new Map());

    return () => {
      lineupRef.current?.destroy();
    };
  }, [collections, mergedData, setFilteredRef, setSelectedRef, getRankingBuilders, adjustRankings]);

  React.useEffect(() => {
    lineupRef.current?.data.getRankings().forEach((ranking) => {
      const structureColumn = ranking.flatColumns.find((col) => col instanceof StructureImageColumn) as StructureImageColumn | undefined;
      structureColumn?.setAlign(align || null);
    })
  }, [align]);

  React.useEffect(() => {
    if (lineupRef.current && rankingRef.current && indexMapRef.current) {
      disableLineUpSelectionListener.current = true;
      if (!selected) {
        lineupRef.current.setSelection([]);
      } else if (!isEqual(selected, lineupRef.current.getSelection())) {
        const selectedIndices = Object.values(selected)
          .flat()
          .map((p) => indexMapRef.current!.get(p)!);
        lineupRef.current.setSelection(selectedIndices);
      }
      onSelectionSetRef.current?.(lineupRef.current, rankingRef.current);
      disableLineUpSelectionListener.current = false;
    }
  }, [selected, onSelectionSetRef]);

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
};
