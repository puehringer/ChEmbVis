import * as React from "react";
import GridLayout, { WidthProvider, Layout } from "react-grid-layout";

const SizedGridLayout = WidthProvider(GridLayout);

declare type GridChild = React.ReactElement<{
  key: Layout["i"];
  gridOptions?: Partial<Layout>;
}>;

export const Grid = ({ children: _children }: { children: GridChild | GridChild[] }) => {
  const [layout, setLayout] = React.useState<Layout[] | null>();

  const children = React.useMemo(() => (!_children || Array.isArray(_children) ? _children.flat() : [_children]).filter((c) => c?.key), [
    _children,
  ]);

  React.useEffect(() => {
    const validChildren = children.filter((c) => c?.key != null);
    const missingLayout = validChildren.filter((c) => !layout?.find((l) => l.i === c.key));

    if (missingLayout.length > 0) {
      setLayout([
        ...(layout?.filter((l) => validChildren.find((c) => c.key === l.i)) || []),
        ...missingLayout.map((c, i) => ({
          // @ts-ignore
          i: c.key as string,
          h: 15,
          w: 6,
          x: i % 2 === 0 ? 0 : 6,
          y: 1,
          ...(c?.props?.gridOptions || {}),
        })),
      ]);
    }
  }, [children, layout]);

  return layout ? (
    <SizedGridLayout
      draggableHandle=".react-grid-item-drag-handle"
      cols={12}
      rowHeight={10}
      className="flex-fill"
      // preventCollision={true}
      verticalCompact={true}
      layout={layout}
      compactType="horizontal"
      onLayoutChange={(layout) => {
        window.dispatchEvent(new Event("resize"));
        setLayout(layout);
      }}
    >
      {children.map((c) => (
        <div key={c.key}>{c}</div>
      ))}
    </SizedGridLayout>
  ) : null;
};
